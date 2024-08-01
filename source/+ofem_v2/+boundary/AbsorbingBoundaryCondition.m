classdef AbsorbingBoundaryCondition < handle & ofem_v2.boundary.FixedBoundary
	%AbsorbingBoundaryCondition Summary of this class goes here
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
        normalVector;
        feBd;
        dim;

        M;
	end

    
    methods
        function objtest = AbsorbingBoundaryCondition(feBd,phys,mesh,name,n)
            objtest.normalVector = n;
            objtest.feBd = ofem_v2.elements.loadFE('HCurl_2D_Order_0');
            objtest.dim = objtest.feBd.dim;
			if ischar(name)
				n = size(mesh.bd, 2);
				for i = 1:n
					if isequal(mesh.bd{1,i},name)
						objtest.boundary = unique(sort(mesh.bd{2,i},2),'rows');
						objtest.index = i;
						break;
					end
				end
			else
				objtest.boundary = unique(sort(mesh.bd{2,name},2),'rows');
				objtest.index = name;
			end
			switch mesh.dim
			case 2
				[~,objtest.edges] = ismember(sort(objtest.boundary,2),mesh.ed,'rows');
				objtest.edges = unique(objtest.edges);
				objtest.nodes = unique(objtest.boundary(:));
			case 3
				[~,objtest.faces] = ismember(sort(objtest.boundary,2),mesh.fa,'rows');
				objtest.faces = unique(objtest.faces);
				objtest.edges = mesh.fa2ed(objtest.faces,:);
				objtest.edges = unique(objtest.edges(:));
				objtest.nodes = unique(objtest.boundary(:));
            end
        end


        function M = assembleMass(objtest, physicalProblem,func)
			[w,l] = ofem_v2.tools.gaussSimplex(objtest.dim,objtest.feBd.degreeMass);
			
            [~,idx] = ismember(objtest.faces,physicalProblem.geometry.el2fa(:));
			elIdx = mod(idx-1,physicalProblem.geometry.Nint)+1;
			faIdx = floor(idx/physicalProblem.geometry.Nint)+1;
			refTet = physicalProblem.geometry.refTet(elIdx);
			refTet = refTet.*faIdx;
            DinvT = physicalProblem.geometry.DinvT(:,:,refTet);

            %calculate the edges on boundary for basis function evaluation
            %very annoying 
            edgeIdx = physicalProblem.geometry.fa2ed(objtest.faces,:);
            edgeIdx = repmat(edgeIdx(:), [1,size(physicalProblem.geometry.el2ed,2)]);
            % nodesofEdge = physicalProblem.geometry.ed(edgeIdx(:),:);
            % nodesofEdge = [physicalProblem.geometry.ed(edgeIdx(:,1),:) physicalProblem.geometry.ed(edgeIdx(:,2),:) physicalProblem.geometry.ed(edgeIdx(:,3),:)];%reshape(nodesofEdge,[length(objtest.faces),6])
            [isEdge,~]=ismember(physicalProblem.geometry.el2ed(elIdx,:),edgeIdx);
            isEdge = repmat(reshape(isEdge',[1,size(physicalProblem.geometry.el2ed,2), length(objtest.faces)]),[3,1]);

			N = max(physicalProblem.DOFs.DOFs);
			dofs = [];
			if ~isempty(physicalProblem.DOFs.n2DOF)
				dofs = [dofs;physicalProblem.geometry.fa(objtest.faces,:)];
			end
			if ~isempty(physicalProblem.DOFs.e2DOF)
				dofs = [dofs,physicalProblem.DOFs.e2DOF(physicalProblem.geometry.fa2ed(objtest.faces,:))];
            end
			if ~isempty(physicalProblem.DOFs.f2DOF)
				dofs = [dofs;physicalProblem.DOFs.f2DOF(objtest.faces)];
			end
			
			Nl     = size(l  ,2); % number of barycentric coordinates
			Nq     = size(w  ,1);
			Ns     = objtest.feBd.DOFsPerElement;
			Nf     = size(physicalProblem.geometry.el,1);
            Ne     = length(refTet);

            v1 = physicalProblem.geometry.co(:,:,physicalProblem.geometry.fa(objtest.faces,1));
            v2 = physicalProblem.geometry.co(:,:,physicalProblem.geometry.fa(objtest.faces,2));
            v3 = physicalProblem.geometry.co(:,:,physicalProblem.geometry.fa(objtest.faces,3));

            vec1 = v2-v1;
            vec2 = v3-v1;
            vec1(3,:,:) = [];
            vec2(3,:,:) = [];

            Dk = [vec1,vec2];

            detD = dot(rot(vec1),vec2,1);

            DinvT = [ -rot(vec2), rot(vec1) ]./detD;

			% faceco gives global quadrature points => eval there
			%faceco = reshape(physicalProblem.geometry.co(:,:,obj.bd(:,1:Nl)'),[],Nl,Nf);
			
			%M      = ofem_v2.tools.matrixarray(zeros(Ns,Ns,Nf));
			M=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));

			if isa(objtest.value,'function_handle')
				% for q  = 1:Nq
				% 	cnt = ones(size(l(:,q),1),1);
				% 	lTemp = mat2cell(l(:,q),cnt);
				% 	phi = objtest.feBd.phi{1}(lTemp{:});
				% 	val = objtest.value(objtest.normalVector);
				% 	F = F + val'*(w(q)*phi);
				% end
			else
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
                    phi(:,:,1) = objtest.feBd.N{1}(lTemp{:});
                    phi(:,:,2) = objtest.feBd.N{2}(lTemp{:});
                    phi =  ofem_v2.tools.matrixarray(phi(:,:,refTet));
                    phi = circshift(phi,1,1);
                    phi = DinvT*phi;
                    % phi(:,phi(1,:,1)==0,:) = [];
                    % normal_reshaped = repmat(reshape(objtest.normalVector, [3, 1, 1]), [1, size(phi,2), size(phi,3)]);
                    % M = M+w(q)*pagemtimes(double(phi),'transpose',double(func*phi),'none');
					M = M + w(q)*pagemtimes(double(func*phi),'transpose',double(phi),'none');
                    clear phi;
				end
            end
            dofs = physicalProblem.geometry.fa2ed(objtest.faces,:);
			M = M*abs(detD);
			I = repmat(dofs,1,Ns)';
            I = I(:);
			J = repelem(dofs,1,Ns)';
            J = J(:);
			objtest.M = sparse(I(:),J(:),M(:),N,N);
			M = objtest.M;
		end
            % 			figure
            %             spy(M)
        end
end    