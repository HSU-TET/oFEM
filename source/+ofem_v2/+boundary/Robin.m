classdef Robin < handle & ofem_v2.boundary.MixedBoundary
	%ROBIN Robin boundary condition  k du/dn + alpha*u = beta*g  for H1
	%elements
	%
	%   bc = Robin(value, alpha, beta, feBd, boundaryName, mesh)
	%
	%   Adds int alpha*u*v ds to M_robin and int beta*g*v ds to the right
	%   hand side. value (g), alpha and beta are constants or function
	%   handles f(X) / f(X,n) of the points X (dim x 1 x N) and the outer
	%   unit normals n (dim x 1 x N), returning 1 x 1 x N.
	%
	%   The trace of the volume basis is integrated on every boundary cell
	%   (traceData / traceBasis), so any degree works in 2D and 3D; feBd is
	%   no longer needed and only kept for compatibility.
	%
	%   lumped = true gives the row-sum lumped boundary mass matrix, which
	%   is what the previous implementation computed (its row and column
	%   indices coincided). The default is the consistent mass matrix.

	properties
		value;
		boundary;
		index;
		name;
		meas;
		normalVector;
		dim;
		feBd;
		edges;
		nodes;
		faces;
		interiors = [];
		elems;
		alpha;
		beta;
		lumped = false;
		M;
		b;
	end

	methods
		function obj = Robin(value, alpha, beta, feBd, boundaryName, mesh)
			obj.value = value;
			obj.alpha = alpha;
			obj.beta = beta;
			obj.feBd = feBd;
			if ischar(boundaryName) || isstring(boundaryName)
				obj.index = find(cellfun(@(s) isequal(s,char(boundaryName)), mesh.bd(1,:)),1);
				if isempty(obj.index)
					error('ofem_v2:boundary:Robin','Boundary ''%s'' not found.',boundaryName);
				end
			else
				obj.index = boundaryName;
			end
			obj.name = mesh.bd{1,obj.index};
			obj.boundary = unique(sort(mesh.bd{2,obj.index},2),'rows');
			obj.nodes = unique(obj.boundary(:));
			obj.elems = size(obj.boundary,1);
			obj.dim = mesh.dim - 1;
			switch mesh.dim
				case 2
					[~,obj.edges] = ismember(obj.boundary,mesh.ed,'rows');
					obj.faces = [];
				case 3
					[~,obj.faces] = ismember(obj.boundary,mesh.fa,'rows');
					obj.edges = unique(reshape(mesh.fa2ed(obj.faces,:),[],1));
			end
			td = ofem_v2.boundary.traceData(mesh, obj.boundary, 1);
			obj.meas = td.meas;
			obj.normalVector = td.normal;
		end

		function b = loadVector(obj, phys)
			fe = phys.element;
			anyFun = isa(obj.value,'function_handle') || isa(obj.beta,'function_handle');
			td = ofem_v2.boundary.traceData(phys.geometry, obj.boundary, ...
				ofem_v2.boundary.Robin.order(fe.degree, anyFun, false));
			P = ofem_v2.boundary.traceBasis(td, phys.geometry, fe.phi);     % 1 x Ns x Nq x Nf
			[~,Ns,Nq,Nf] = size(P);
			g = obj.evaluate(obj.value, td).*obj.evaluate(obj.beta, td);   % 1 x 1 x Nq x Nf
			jw = reshape(td.w,1,1,Nq).*reshape(td.jac,1,1,1,Nf);
			F = reshape(sum(P.*g.*jw,3),Ns,Nf);
			dofs = phys.DOFs.el2DOF(td.elem,:);
			obj.b = sparse(reshape(dofs',[],1),1,F(:),phys.DOFs.Nd,1);
			b = obj.b;
		end

		function M = assembleMass(obj, phys)
			fe = phys.element;
			td = ofem_v2.boundary.traceData(phys.geometry, obj.boundary, ...
				ofem_v2.boundary.Robin.order(fe.degree, isa(obj.alpha,'function_handle'), true));
			P = ofem_v2.boundary.traceBasis(td, phys.geometry, fe.phi);     % 1 x Ns x Nq x Nf
			[~,Ns,Nq,Nf] = size(P);
			a = obj.evaluate(obj.alpha, td);
			jw = reshape(td.w,1,1,Nq).*reshape(td.jac,1,1,1,Nf).*a;
			P = reshape(P,Ns,Nq,Nf);
			Me = pagemtimes(P.*reshape(jw,1,Nq,Nf),'none',P,'transpose');  % Ns x Ns x Nf
			dofs = phys.DOFs.el2DOF(td.elem,:);
			I = repmat(dofs,1,Ns)';
			J = repelem(dofs',Ns,1);
			N = phys.DOFs.Nd;
			M = sparse(I(:),J(:),Me(:),N,N);
			if obj.lumped
				M = spdiags(full(sum(M,2)),0,N,N);
			end
			obj.M = M;
		end
	end

	methods (Access = private)
		function v = evaluate(~, f, td)
			% values at the quadrature points: 1 x 1 x Nq x Nf
			[~,Nq,Nf] = size(td.X);
			if ~isa(f,'function_handle')
				v = f*ones(1,1,Nq,Nf);
				return
			end
			v = zeros(1,1,Nq,Nf);
			for q = 1:Nq
				if nargin(f) >= 2
					fq = f(td.X(:,q,:),td.normal);
				else
					fq = f(td.X(:,q,:));
				end
				v(1,1,q,:) = reshape(double(fq),1,1,1,Nf);
			end
		end
	end

	methods (Static)
		function q = order(p, anyFun, mass)
			% exact for constant data; for function data up to degree p+2
			if anyFun
				q = 2*p + 2;
			elseif mass
				q = 2*p;
			else
				q = p;
			end
		end
	end
end
