classdef WaveguideInput < handle & ofem_v2.boundary.FixedBoundary
	%WAVEGUIDEINPUT Summary of this class goes here
	%   Detailed explanation goes here
	
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
			% handle edges
			% TODO: 1 dim quadrature
			[w,l] = ofem_v2.tools.gaussSimplex(1,1);
			eDofs = phys.DOFs.e2DOF(obj.edges,:);
			eDofs = reshape(eDofs',1,phys.element.degree+1,[]);
			v1 = phys.geometry.co(:,:,phys.geometry.ed(obj.edges,1));
			v2 = phys.geometry.co(:,:,phys.geometry.ed(obj.edges,2));
			tang = v2-v1;
			tang = double(tang);%*ofem_v2.tools.matrixarray(1/vecnorm(tang));

			% Interpolate at edge center using stefans function
			obj.cell2edge(phys);

			b = ofem_v2.tools.matrixarray(zeros(size(eDofs)));
			%obj.b = ofem_v2.tools.matrixarray(reshape(obj.b,1,size(eDofs,2),[]));

			for q = 1:length(w)
				b = b + pagemtimes(obj.value,'transpose',tang,'none');
			end

			obj.u = sparse(eDofs(:),1,b(:),phys.DOFs.Nd,1);

			u = obj.u;
			phys.DOFs.reduceDOFs(obj);
		end

		function cell2edge(obj, phys)
			% convert cell data from modes to edge data
			v1 = phys.geometry.co(:,:,phys.geometry.ed(obj.edges,1));
			v2 = phys.geometry.co(:,:,phys.geometry.ed(obj.edges,2));
			xq = (v1+v2)/2;
			xq = xq(1:2,:,:);
			n = 10;%5;		
    		co_idx = knnsearch(obj.x',squeeze(xq)','K',n)';
		
    		ui_x = (sum((reshape(obj.u_in(1,co_idx(:,:)'),[],n)),2)/n);
    		ui_y = (sum((reshape(obj.u_in(2,co_idx(:,:)'),[],n)),2)/n);
			obj.value = reshape([ui_x,ui_y]',2,1,[]);
			obj.value(3,:,:) = 0;
		end
	end
end

