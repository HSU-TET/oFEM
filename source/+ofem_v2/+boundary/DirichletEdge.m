classdef DirichletEdge < handle & ofem_v2.boundary.FixedBoundary
	%DIRICHLETEDGE Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		value;
		boundary;
        nodes;
		edges;
		faces;
        interiors;
		index;
		boundaryName;
		
		u;
	end
	
	methods
		function obj = DirichletEdge(value, boundaryName, mesh)
			if isscalar(value)
				error('ofem_v2:boundaries:DirichletEdge',...
						'Dirichlet boundary needs to be a vector');
			end
            obj.value = value;
            
			if ischar(boundaryName)
				n = size(mesh.bd, 2);
				for i = 1:n
					if isequal(mesh.bd{1,i},boundaryName)
						obj.boundary = unique(sort(mesh.bd{2,i},2),'rows');
						obj.index = i; 
						break;
					end
				end
			else
				obj.boundary = unique(sort(mesh.bd{2,boundaryName},2),'rows');
				obj.index = boundaryName;
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
            % handle edges
            % TODO: 1 dim quadrature
            [w,l] = ofem_v2.tools.gaussSimplex(1,phys.element.degreeStiff/2);
            eDofs = phys.DOFs.e2DOF(obj.edges,:);
            eDofs = reshape(eDofs',1,phys.element.degree+1,[]);
            v1 = phys.geometry.co(:,:,phys.geometry.ed(obj.edges,1));
            v2 = phys.geometry.co(:,:,phys.geometry.ed(obj.edges,2));
            tang = v2-v1;
            tang = tang*ofem_v2.tools.matrixarray(1/vecnorm(tang));
            
            b = ofem_v2.tools.matrixarray(zeros(size(eDofs)));
            %obj.b = ofem_v2.tools.matrixarray(reshape(obj.b,1,size(eDofs,2),[]));
            
            for q = 1:length(w)
				if phys.element.dim == 3
					phi = phys.element.N{1}(l(q),0,0);
				else
					phi = phys.element.N{1}(l(q),0);
				end
                phi = phi(1,1:phys.element.degree+1);
                b = b + (obj.value'*tang)*phi;
            end
            
            obj.u = sparse(eDofs(:),1,b(:),phys.DOFs.Nd,1);
            
            % handle faces
            % TODO: 2 dim quadrature
            if phys.element.degree > 1 && phys.element.dim == 3
                [w,l] = ofem_v2.tools.gaussSimplex(2,phys.element.degreeStiff/2);

                fDofs = phys.DOFs.f2DOF(obj.faces,:);
                fDofs = reshape(fDofs',1,size(fDofs,2),[]);
                
                eoff = phys.element.edgeDOFs;
                foff = size(fDofs,2);

                v1 = phys.geometry.co(:,:,phys.geometry.fa(obj.faces,1));
                v2 = phys.geometry.co(:,:,phys.geometry.fa(obj.faces,2)); 
                v3 = phys.geometry.co(:,:,phys.geometry.fa(obj.faces,3));

                A = vecnorm(cross(v2-v1,v3-v1))/2;
                
                b = ofem_v2.tools.matrixarray(zeros(size(fDofs)));
                
                % Find the elements containing the faces (annoying)
                
                [~,idx] = ismember(obj.faces,phys.geometry.el2fa(:));
                elIdx = mod(idx-1,phys.geometry.Nint)+1;
                faIdx = floor(idx/phys.geometry.Nint)+1;
                Dk = phys.geometry.Dk(:,:,elIdx);
                refTet = phys.geometry.refTet(elIdx);
                refTet = refTet.*faIdx;
                soff = [eoff+1:eoff+foff;eoff+1+foff:eoff+2*foff;eoff+1+2*foff:eoff+3*foff;eoff+1+3*foff:eoff+4*foff];
                for q = 1:length(w)
                    phi = phys.element.N{1}(l(1,q),l(2,q),0);
                    phitot = phi(:,soff(1,:));
                    phi = phys.element.N{2}(l(1,q),l(2,q),0);
                    phitot(:,:,2) = phi(:,soff(1,:));
                    phi = phys.element.N{1}(l(1,q),0,l(2,q));
                    phitot(:,:,3) = phi(:,soff(2,:));
                    phi = phys.element.N{2}(l(1,q),0,l(2,q));
                    phitot(:,:,4) = phi(:,soff(2,:));
                    phi = phys.element.N{1}(0,l(1,q),l(2,q));
                    phitot(:,:,5) = phi(:,soff(3,:));
                    phi = phys.element.N{2}(0,l(1,q),l(2,q));
                    phitot(:,:,6) = phi(:,soff(3,:));
                    phi = phys.element.N{1}(l(1,q),l(2,q),1-l(1,q)-l(2,q));
                    phitot(:,:,7) = phi(:,soff(4,:));
                    phi = phys.element.N{2}(l(1,q),l(2,q),1-l(1,q)-l(2,q));
                    phitot(:,:,8) = phi(:,soff(4,:));
                    
                    phi = Dk*ofem_v2.tools.matrixarray(phitot(:,:,refTet));
                    
                    b = b + obj.value'*phi*ofem_v2.tools.matrixarray(1/A);
                    
                end
                obj.u = obj.u + sparse(fDofs(:),1,b(:),phys.DOFs.Nd,1);
            end
            
            u = obj.u;
			phys.DOFs.reduceDOFs(obj);
        end
	end
end

























