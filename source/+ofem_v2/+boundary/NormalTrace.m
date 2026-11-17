classdef NormalTrace < handle & ofem_v2.boundary.FixedBoundary
	%NORMALTRACE Homogeneous essential boundary condition u.n = 0 for
	%H(div) elements.
	%
	%   bc = NormalTrace(name, mesh) fixes all DOFs on the boundary faces
	%   (dim = 3) or boundary edges (dim = 2) of the boundary 'name' (or
	%   its index). Use with ofem_v2.elements.HDivElement.

	properties
		value = 0;
		boundary;
		nodes = [];
		edges = [];
		faces = [];
		interiors = [];
		index;
		name;
		u;
	end

	methods
		function obj = NormalTrace(name, mesh)
			if ischar(name) || isstring(name)
				obj.name = char(name);
				obj.index = find(cellfun(@(s) isequal(s,obj.name), mesh.bd(1,:)),1);
				if isempty(obj.index)
					error('ofem_v2:boundary:NormalTrace', ...
						'Boundary ''%s'' not found.', obj.name);
				end
			else
				obj.index = name;
				obj.name = mesh.bd{1,name};
			end
			obj.boundary = unique(sort(mesh.bd{2,obj.index},2),'rows');
			switch mesh.dim
				case 2
					[~,obj.edges] = ismember(obj.boundary,mesh.ed,'rows');
					obj.edges = unique(obj.edges);
				case 3
					[~,obj.faces] = ismember(obj.boundary,mesh.fa,'rows');
					obj.faces = unique(obj.faces);
			end
		end

		function u = loadVector(obj, phys)
			% the prescribed normal trace is zero
			obj.u = sparse(phys.DOFs.Nd,1);
			u = obj.u;
			phys.DOFs.reduceDOFs(obj);
		end
	end
end
