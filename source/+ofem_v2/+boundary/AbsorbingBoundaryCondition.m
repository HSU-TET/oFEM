classdef AbsorbingBoundaryCondition < handle & ofem_v2.boundary.MixedBoundary
	%ABSORBINGBOUNDARYCONDITION first order absorbing boundary for H(curl)
	%problems
	%
	%   bc = AbsorbingBoundaryCondition(feBd, phys, mesh, name, n, factor)
	%
	%   Adds factor * int (n x u).(n x v) ds (the tangential boundary mass
	%   matrix, e.g. factor = 1i*k for the Silver-Mueller condition) to
	%   M_robin; there is no right hand side. factor is a (complex) constant
	%   or a function handle f(X) / f(X,n) returning 1 x 1 x N.
	%
	%   The trace of the volume basis is integrated on every boundary cell
	%   (traceData / traceBasis), so any degree works in 2D and 3D and the
	%   faces may have any orientation. feBd, phys and n are not needed
	%   any more and only kept for compatibility (the outer normals are
	%   computed from the mesh). The previous version handled lowest order
	%   elements on faces in the xy plane only.

	properties
		value = 0;
		boundary;
		nodes = [];
		edges;
		faces;
		interiors = [];
		index;
		name;
		normalVector;
		feBd;
		dim;
		factor;
		M;
	end

	methods
		function obj = AbsorbingBoundaryCondition(feBd, phys, mesh, name, n, factor) %#ok<INUSL>
			obj.feBd = feBd;
			obj.factor = factor;
			if ischar(name) || isstring(name)
				obj.index = find(cellfun(@(s) isequal(s,char(name)), mesh.bd(1,:)),1);
				if isempty(obj.index)
					error('ofem_v2:boundary:AbsorbingBoundaryCondition','Boundary ''%s'' not found.',name);
				end
			else
				obj.index = name;
			end
			obj.name = mesh.bd{1,obj.index};
			obj.boundary = unique(sort(mesh.bd{2,obj.index},2),'rows');
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
			obj.normalVector = td.normal;
		end

		function M = assembleMass(obj, phys)
			fe = phys.element;
			isFun = isa(obj.factor,'function_handle');
			td = ofem_v2.boundary.traceData(phys.geometry, obj.boundary, 2*max(fe.degree,1) + 2*isFun);
			B = ofem_v2.boundary.traceBasis(td, phys.geometry, fe.N, 'covariant');  % dim x Ns x Nq x Nf
			[dim,Ns,Nq,Nf] = size(B);
			n = reshape(td.normal,dim,1,1,Nf);
			if dim == 2
				T = sum([-n(2,:,:,:); n(1,:,:,:)].*B,1);
			else
				T = B - n.*sum(n.*B,1);
			end
			nt = size(T,1);
			jw = reshape(td.w,1,1,Nq).*reshape(td.jac,1,1,1,Nf);
			if isFun
				for q = 1:Nq
					if nargin(obj.factor) >= 2
						fq = obj.factor(td.X(:,q,:),td.normal);
					else
						fq = obj.factor(td.X(:,q,:));
					end
					jw(1,1,q,:) = jw(1,1,q,:).*reshape(fq,1,1,1,Nf);
				end
				fac = 1;
			else
				fac = obj.factor;
			end
			Tm = reshape(permute(T,[1,3,2,4]),nt*Nq,Ns,Nf);
			Tw = reshape(permute(T.*jw,[1,3,2,4]),nt*Nq,Ns,Nf);
			Me = fac*pagemtimes(Tm,'transpose',Tw,'none');
			dofs = phys.DOFs.el2DOF(td.elem,:);
			I = repmat(dofs,1,Ns)';
			J = repelem(dofs',Ns,1);
			obj.M = sparse(I(:),J(:),Me(:),phys.DOFs.Nd,phys.DOFs.Nd);
			M = obj.M;
		end

		function b = loadVector(obj, phys) %#ok<INUSL>
			b = sparse(phys.DOFs.Nd,1);
		end
	end
end
