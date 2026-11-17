classdef DirichletEdge < handle & ofem_v2.boundary.FixedBoundary
	%DIRICHLETEDGE prescribed tangential trace n x u = n x g of an H(curl)
	%field on a boundary
	%
	%   bc = DirichletEdge(value, name, mesh)
	%
	%   value is a constant vector, 0, or a function handle g(X) / g(X,n).
	%   The trace is set by an L2 projection onto the boundary DOFs (see
	%   tangentialProjection), the DOFs of the boundary edges and faces are
	%   fixed.

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

	methods
		function obj = DirichletEdge(value, name, mesh)
			% value: constant vector (dim x 1), zero, or a function handle
			% g(X) / g(X,n) of the points X (dim x 1 x N) returning dim x 1 x N
			if ~isa(value,'function_handle') && isscalar(value) && value ~= 0
				error('ofem_v2:boundaries:DirichletEdge',...
					'Dirichlet boundary needs to be a vector or a function handle');
			end
			obj.value = value;

			if ischar(name)
				n = size(mesh.bd, 2);
				for i = 1:n
					if isequal(mesh.bd{1,i},name)
						obj.boundary = unique(sort(mesh.bd{2,i},2),'rows');
						obj.index = i;
						break;
					end
				end
			else
				obj.boundary = unique(sort(mesh.bd{2,name},2),'rows');
				obj.index = name;
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

		function u = loadVector(obj, phys)
			% L2 projection of the tangential trace onto the DOFs of the
			% boundary edges (and faces in 3D): exact for data in the
			% discrete trace space, any degree. The DOFs are fixed.
			% (The previous version summed value.t once per quadrature
			% point without weights, basis or edge length, which is only
			% correct for zero data.)
			obj.u = ofem_v2.boundary.tangentialProjection(phys, obj.boundary, ...
				obj.edges, obj.faces, obj.value);
			u = obj.u;
			phys.DOFs.reduceDOFs(obj);
		end
	end
end

























