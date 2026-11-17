classdef WaveguideInput < handle & ofem_v2.boundary.FixedBoundary
	%WAVEGUIDEINPUT tangential field of a waveguide mode on a port (3D)
	%
	%   bc = WaveguideInput(u_in, x, name, mesh)
	%
	%   u_in  transverse field samples (2 x Np) at the points x (2 x Np) of
	%         the port plane (coordinates 1 and 2 of the mesh points)
	
	properties
		value;
		boundary;
		nodes;
		edges;
		faces;
		interiors;
		index;
		name;

		u_in;
		x;

		u;
	end
	
	methods
		function obj = WaveguideInput(u_in, x, name, mesh)
% 			if size(u_in,1) ~= 2 || size(u_in,2) ~= 2
% 				error('ofem_v2:boundaries:WaveguideInput',...
% 					'Waveguide Input has to be an array of vectors size(2,numEdges)');
% 			end
			obj.u_in = u_in;
			obj.x = x;

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
					error('ofem_v2:boundaries:WaveguideInput',...
					'Waveguide Input is only defined in 3D');
				case 3
					[~,obj.faces] = ismember(sort(obj.boundary,2),mesh.fa,'rows');
					obj.faces = unique(obj.faces);
					obj.edges = mesh.fa2ed(obj.faces,:);
					obj.edges = unique(obj.edges(:));
					obj.nodes = unique(obj.boundary(:));
			end
		end

		function u = loadVector(obj, phys)
			% The mode field at a point is the mean of the n nearest samples
			% u_in at x (coordinates 1:2 of the point, i.e. a port in a plane
			% z = const); it is set by the L2 projection of its tangential
			% trace (see tangentialProjection). The previous version assigned
			% value.t to every edge DOF, which has the wrong sign for the
			% Whitney DOFs and is wrong for the higher order DOFs.
			obj.u = ofem_v2.boundary.tangentialProjection(phys, obj.boundary, ...
				obj.edges, obj.faces, @(X) obj.mode(X));
			u = obj.u;
			phys.DOFs.reduceDOFs(obj);
		end

		function v = mode(obj, X)
			% interpolated mode field at the points X (3 x 1 x N): 3 x 1 x N
			n = min(obj.k, size(obj.x,2));
			xq = reshape(double(X(1:2,1,:)),2,[])';
			idx = knnsearch(obj.x',xq,'K',n);
			ux = mean(reshape(obj.u_in(1,idx),[],n),2);
			uy = mean(reshape(obj.u_in(2,idx),[],n),2);
			v = reshape([ux, uy, zeros(size(ux))]',3,1,[]);
		end
	end
end

