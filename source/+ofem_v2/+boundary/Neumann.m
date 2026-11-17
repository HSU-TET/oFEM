classdef Neumann < handle & ofem_v2.boundary.NaturalBoundary
	%NEUMANN natural boundary condition for H1 elements: adds int g*v ds
	%to the right hand side.
	%
	%   bc = Neumann(value, feBd, boundaryName, mesh)
	%
	%   value         constant g, or a function handle g(X) or g(X,n) of
	%                 the quadrature points X (dim x 1 x Nf) and the outer
	%                 unit normals n (dim x 1 x Nf)
	%   feBd          H1 element on the boundary (dimension dim-1) of the
	%                 same order as the volume element; in 3D only order 1
	%   boundaryName  name (or index) of the boundary in the mesh

	properties
		value;
		boundary;       % boundary edges/faces as sorted node lists
		index;
		meas;           % length/area per boundary entity (1 x 1 x Nf)
		normalVector;   % outer unit normal per boundary entity (dim x 1 x Nf)
		dim;            % dimension of the boundary
		feBd;
		edges;          % global edges of the boundary entities (2D: aligned with boundary)
		nodes;
		faces;          % global faces (3D: aligned with boundary)
		interiors = [];
		elems;          % number of boundary entities
		b;
	end

	methods
		function obj = Neumann(value, feBd, boundaryName, mesh)
			obj.value = value;
			obj.feBd = feBd;
			if ischar(boundaryName) || isstring(boundaryName)
				obj.index = find(cellfun(@(s) isequal(s,char(boundaryName)), mesh.bd(1,:)),1);
				if isempty(obj.index)
					error('ofem_v2:boundary:Neumann','Boundary ''%s'' not found.',boundaryName);
				end
			else
				obj.index = boundaryName;
			end
			obj.boundary = unique(sort(mesh.bd{2,obj.index},2),'rows');
			obj.nodes = unique(obj.boundary(:));
			obj.elems = size(obj.boundary,1);
			x1 = mesh.co(:,:,obj.boundary(:,1));
			x2 = mesh.co(:,:,obj.boundary(:,2));
			switch mesh.dim
				case 2
					obj.dim = 1;
					[~,obj.edges] = ismember(obj.boundary,mesh.ed,'rows');
					t = x2-x1;
					obj.meas = sqrt(sum(t.^2,1));
					n = [t(2,:,:); -t(1,:,:)]./obj.meas;
					adj = adjacentElements(mesh.el2ed,obj.edges);
				case 3
					obj.dim = 2;
					[~,obj.faces] = ismember(obj.boundary,mesh.fa,'rows');
					obj.edges = unique(reshape(mesh.fa2ed(obj.faces,:),[],1));
					x3 = mesh.co(:,:,obj.boundary(:,3));
					n = cross(x2-x1,x3-x1,1);
					obj.meas = sqrt(sum(n.^2,1))/2;
					n = n./(2*obj.meas);
					adj = adjacentElements(mesh.el2fa,obj.faces);
			end
			% orient the normals outwards: away from the adjacent element
			Nv = mesh.dim+1;
			c = mean(reshape(mesh.co(:,:,mesh.el(adj,1:Nv)'),mesh.dim,Nv,[]),2);
			s = sign(sum(n.*(x1-c),1));
			s(s == 0) = 1;
			obj.normalVector = n.*s;
		end

		function b = loadVector(obj, phys)
			fe = obj.feBd;
			p = fe.degree;
			if obj.dim == 2 && p > 1
				error('ofem_v2:boundary:Neumann', ...
					'Neumann boundaries in 3D are implemented for order 1 only.');
			end
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim,p+2);
			if obj.dim == 1
				l = l(:)'; w = w(:);
			end
			g = phys.geometry;
			Nf = obj.elems;
			Nv = obj.dim+1;
			xv = cell(1,Nv);
			for k = 1:Nv
				xv{k} = g.co(:,:,obj.boundary(:,k));
			end
			Ns = fe.DOFsPerElement;
			F = zeros(1,Ns,Nf);
			for q = 1:numel(w)
				c = num2cell(l(:,q));
				phi = fe.phi{1}(c{:});
				if isa(obj.value,'function_handle')
					X = xv{1};
					for k = 2:Nv
						X = X + l(k-1,q)*(xv{k}-xv{1});
					end
					if nargin(obj.value) >= 2
						val = obj.value(X,obj.normalVector);
					else
						val = obj.value(X);
					end
					val = reshape(double(val),1,1,Nf);
				else
					val = obj.value;
				end
				F = F + w(q)*val.*phi;
			end
			% w sums to the reference measure (1 for edges, 1/2 for triangles)
			F = F.*obj.meas*factorial(obj.dim);

			% DOFs: vertex functions, then the edge functions (2D)
			D = phys.DOFs;
			dofs = zeros(Nf,Ns);
			for k = 1:Nv
				dofs(:,k) = D.n2DOF(obj.boundary(:,k));
			end
			if Ns > Nv
				dofs(:,Nv+1:end) = D.e2DOF(obj.edges,:);
			end
			obj.b = sparse(reshape(dofs',[],1),1,F(:),D.Nd,1);
			b = obj.b;
		end
	end
end

function adj = adjacentElements(el2x,idx)
	% one element containing each of the entities idx
	Ne = size(el2x,1);
	owner = zeros(max(el2x(:)),1);
	owner(el2x(:)) = repmat((1:Ne)',size(el2x,2),1);
	adj = owner(idx);
end
