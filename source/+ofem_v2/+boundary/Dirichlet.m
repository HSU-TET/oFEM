classdef Dirichlet < handle
	%   Summary of this class goes here
	%   Detailed explanation goes here
	%% Proteries
	properties
		value;
		boundary;
		nodes;
		edges;
		faces;
		interiors;
		index;
		name;
		
		u;
		
	end
	%% Methodes
	methods
		function obj = Dirichlet(value, boundaryName, mesh)
			obj.value = value;
			
			n = size(mesh.bd, 2);
			for i = 1:n
				if isequal(mesh.bd{1,i},boundaryName)
					obj.boundary = unique(sort(mesh.bd{2,i},2),'rows');
					obj.index = i;
					obj.name = boundaryName;
					break;
				end
			end
			switch mesh.dim
				case 2
					[~,obj.edges] = ismember(sort(obj.boundary,2),mesh.ed,'rows');
					obj.edges = unique(obj.edges);
					obj.nodes = unique(obj.boundary(:));
				case 3
					[~,obj.faces] = ismember(sort(obj.boundary,2),mesh.fa,'rows');
					obj.faces = unique(obj.faces);
					obj.edges = mesh.fa2ed(obj.faces,:);
					obj.edges = unique(obj.edges(:));
					obj.nodes = unique(obj.boundary(:));
			end
		end
		
		function u = loadVector(obj, physicalProblem)
			N = max(physicalProblem.DOFs.DOFs);
			dofs = [];
			if ~isempty(physicalProblem.DOFs.n2DOF)
				dofs = [dofs;physicalProblem.DOFs.n2DOF(obj.nodes)];
			end
			if ~isempty(physicalProblem.DOFs.e2DOF)
				dofs = [dofs;physicalProblem.DOFs.e2DOF(obj.edges)];
			end
			if ~isempty(physicalProblem.DOFs.f2DOF)
				dofs = [dofs;physicalProblem.DOFs.f2DOF(obj.faces)];
			end
			obj.u = sparse(dofs,1,obj.value,N,1);
			
			u = obj.u;
			physicalProblem.DOFs.reduceDOFs(obj);
			
		end
	end
end























