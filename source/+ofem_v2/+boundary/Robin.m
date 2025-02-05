classdef Robin < handle & ofem_v2.boundary.MixedBoundary
	%ROBINBOUNDARY Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		value;
		boundary;
		index;
		meas;
		normalVector;
		bd;
		dim;
		feBd;
		edges;
		nodes;
		faces;
		elems;
		alpha;
		beta;
		
		M;
		b;
	end
	
	methods
		function obj=Robin(value, alpha, beta, feBd, boundaryName, mesh)
			obj.value = value;
			obj.feBd = feBd;
			obj.alpha = alpha;
			obj.beta = beta;
			
			n = size(mesh.bd, 2);
			for i = 1:n
				if isequal(mesh.bd{1,i},boundaryName)
					obj.boundary = mesh.bd{2,i};%unique(sort(mesh.bd{2,i},2),'rows');
					obj.index= i;
					break;
				end
			end
			
			switch mesh.dim
				case 2
					[~,obj.edges] = ismember(sort(obj.boundary,2),mesh.ed,'rows');
					obj.edges = unique(obj.edges);
					obj.nodes = unique(obj.boundary(:));
					edges = mesh.co(:,:,mesh.bd{2, obj.index}(:,2))- mesh.co(:,:,mesh.bd{2, obj.index}(:,1));
					obj.meas = ofem_v2.tools.matrixarray(sqrt(dot(edges,edges,1)));
					obj.normalVector = -edges.rot*(1/obj.meas);
					obj.dim = 1;
					obj.elems = size(obj.edges,1);
				case 3
					[~,obj.faces] = ismember(sort(obj.boundary,2),mesh.fa,'rows');
					obj.faces = unique(obj.faces);
					obj.edges = mesh.fa2ed(obj.faces,:);
					obj.edges = unique(obj.edges(:));
					obj.nodes = unique(obj.boundary(:));
					e1 = mesh.co(:,:,mesh.bd{2, obj.index}(:,2))-mesh.co(:,:,mesh.bd{2, obj.index}(:,1));
					e2 = mesh.co(:,:,mesh.bd{2, obj.index}(:,3))-mesh.co(:,:,mesh.bd{2, obj.index}(:,1));
					obj.normalVector = cross(e1,e2,1);
					obj.meas = sqrt(dot(obj.normalVector,obj.normalVector,1))/2;
					obj.dim = 2;
					obj.elems = size(obj.faces,1);
			end
		end
		
		function b = loadVector(obj, physicalProblem)			
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim, obj.feBd.degree);
			
			N = max(physicalProblem.DOFs.DOFs);
			dofs = [];
			if ~isempty(physicalProblem.DOFs.n2DOF)
				dofs = [dofs;physicalProblem.geometry.fa(obj.faces,:)];
			end
			if ~isempty(physicalProblem.DOFs.e2DOF)
				dofs = [dofs,physicalProblem.DOFs.e2DOF(physicalProblem.geometry.fa2ed(obj.faces,:))];
			end
			if ~isempty(physicalProblem.DOFs.f2DOF)
				dofs = [dofs;physicalProblem.DOFs.f2DOF(obj.faces)];
			end
			
			Nl     = size(l  ,2); % number of barycentric coordinates
			Nq     = size(w  ,1);
			Ns     = obj.feBd.DOFsPerElement;
			Nf     = obj.elems;
			
			% faceco gives global quadrature points => eval there
			%faceco = reshape(physicalProblem.geometry.co(:,:,obj.bd(:,1:Nl)'),[],Nl,Nf);
			
			F      = ofem_v2.tools.matrixarray(zeros(1,Ns,Nf));
			
			if isa(obj.value,'function_handle')
				for q  = 1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi = obj.feBd.phi{1}(lTemp{:});
					val = obj.value(obj.normalVector);
					F = F + val'*(w(q)*phi);
				end
			else
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi = obj.feBd.phi{1}(lTemp{:});
					F = F + obj.value*(w(q)*phi);
				end
			end
			
			F = F*obj.meas*obj.beta;
			I = dofs';
			obj.b = sparse(I(:),1,F(:),N,1);
			b = obj.b;
		end
		
		function M = assembleMass(obj, physicalProblem)
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim, obj.feBd.degree);
			
			N = max(physicalProblem.DOFs.DOFs);
			dofs = [];
			if ~isempty(physicalProblem.DOFs.n2DOF)
				dofs = [dofs;physicalProblem.geometry.fa(obj.faces,:)];
			end
			if ~isempty(physicalProblem.DOFs.e2DOF)
				dofs = [dofs,physicalProblem.DOFs.e2DOF(physicalProblem.geometry.fa2ed(obj.faces,:))];
            end
			if ~isempty(physicalProblem.DOFs.f2DOF)
				dofs = [dofs;physicalProblem.DOFs.f2DOF(obj.faces)];
			end
			
			Nl     = size(l  ,2); % number of barycentric coordinates
			Nq     = size(w  ,1);
			Ns     = obj.feBd.DOFsPerElement;
			Nf     = obj.elems;
			
			% faceco gives global quadrature points => eval there
			%faceco = reshape(physicalProblem.geometry.co(:,:,obj.bd(:,1:Nl)'),[],Nl,Nf);
			
			M      = ofem_v2.tools.matrixarray(zeros(Ns,Ns,Nf));
			
			if isa(obj.value,'function_handle')
				for q  = 1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi = obj.feBd.phi{1}(lTemp{:});
					val = obj.value(obj.normalVector);
					F = F + val'*(w(q)*phi);
				end
			else
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi = obj.feBd.phi{1}(lTemp{:});
					M = M + obj.alpha*w(q)*(phi'*phi);
				end
			end
			
			I = repmat(dofs,1,Ns)';
            I = I(:);
			J = repelem(dofs,1,Ns)';
            J = J(:);
			M = M*obj.meas;
			obj.M = sparse(I(:),J(:),M(:),N,N);
			M = obj.M;
		end
	end
end






















