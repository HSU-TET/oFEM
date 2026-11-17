classdef Physical_Problem < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        element;
        geometry;
        
        has_stiffness;
        paraS;
        S;
        SPart = {};
        
        has_damping;
        D;
        DPart = {};
        
        has_mass;
        paraM;
        M;
        MPart = {};
        
        F = {};
        b;
		b_r = {};
        
        M_robin;
		M_robinPart;
        
        bc;
        ic;
        
        DOFs;
        
        u;
        
        matlib;

        % elements per call of the element assembly routines; the triplets
        % and the working memory of sparse() exist only for one chunk
        assemblyChunk = 1e6;
    end
    
    methods
        
        
        function obj = Physical_Problem(element, geometry, has_s, has_m, has_d)
            obj.element = element;
            obj.geometry = geometry;
            obj.has_stiffness = has_s;
            obj.has_damping = has_d;
            obj.has_mass = has_m;
            obj.u = zeros(obj.geometry.Nco,1);
            
        end
        
        function setParaS(obj, namePara)
            obj.paraS = namePara;
        end
        
        function setParaM(obj,namePara)
            obj.paraM = namePara;
        end
        
        function setBoundaryCondition(obj, BC)
            obj.geometry.bd{3,BC.index} = BC;
        end
        
        function assemble(obj,quiet)
            %Nc = size(obj.geometry.co  ,3);
			if exist('quiet','var') == 1
				obj.checkBoundaries(1);
			else
				obj.checkBoundaries(2);
			end
            Np = size(obj.geometry.parts, 2);
            
            %obj.DOFs = 1:Nc;
            Nc = obj.DOFs.Nd;
            
            obj.S=sparse(Nc,Nc);
            obj.D=sparse(Nc,Nc);
            obj.M=sparse(Nc,Nc);
            obj.b=sparse(Nc,1);
			
			obj.M_robin = sparse(Nc,Nc);
            
            obj.u = zeros(Nc,1);
            
            obj.M_robin=sparse(Nc,Nc);
            
            has_f = size(obj.geometry.parts,1)>=4;
            
            % start assembling
            tic;
            
            for i=1:Np
                pIdx = obj.geometry.parts{3,i};
                % with a single part S and SPart{1} share their memory
                if obj.has_stiffness
                    obj.SPart{i} = obj.assembleChunked(@(p,m) obj.element.assembleStiffness(obj,p,m), ...
                        pIdx,obj.geometry.parts{2,i}.(obj.paraS));
                    if i == 1, obj.S = obj.SPart{i}; else, obj.S = obj.S + obj.SPart{i}; end
                end
                if obj.has_mass
                    obj.MPart{i} = obj.assembleChunked(@(p,m) obj.element.assembleMass(obj,p,m), ...
                        pIdx,obj.geometry.parts{2,i}.(obj.paraM));
                    if i == 1, obj.M = obj.MPart{i}; else, obj.M = obj.M + obj.MPart{i}; end
                end
                if has_f
                    if ~isempty(obj.geometry.parts{4,i})
                        obj.F{i} = obj.geometry.parts{4,i}.force(obj,pIdx);
                        obj.b = obj.b+obj.F{i};
                    end
                end                
            end
            
            if size(obj.geometry.bd,1) == 3
                nB= size(obj.geometry.bd,2);
                for i=1:nB
                    if ~isempty(obj.geometry.bd{3,i})
						if isa(obj.geometry.bd{3,i},'ofem_v2.boundary.FixedBoundary')
                        	obj.u = obj.u + obj.geometry.bd{3,i}.loadVector(obj);
						elseif isa(obj.geometry.bd{3,i},'ofem_v2.boundary.NaturalBoundary')
							obj.b = obj.b + obj.geometry.bd{3,i}.loadVector(obj);
						elseif isa(obj.geometry.bd{3,i},'ofem_v2.boundary.MixedBoundary')
							obj.b_r{i} = obj.geometry.bd{3,i}.loadVector(obj);
							obj.b = obj.b + obj.b_r{i};
							obj.M_robinPart{i} = obj.geometry.bd{3,i}.assembleMass(obj);
							obj.M_robin = obj.M_robin + obj.M_robinPart{i};
						end
                    end
                end
                if any(obj.u)
                    % lift the Dirichlet values (without forming S+M+D+M_robin)
                    obj.b = obj.b - obj.S*obj.u - obj.M*obj.u - obj.D*obj.u - obj.M_robin*obj.u;
                end
            end
        end

        function A = assembleChunked(obj,fun,pIdx,mat)
            % A = fun(pIdx,mat) assembled in chunks of obj.assemblyChunk
            % elements; per element materials (nc x nc x Ne) are split too.
            % (Column blocks instead of sums are slower: the DOF numbering
            % is not local, elements would be assembled ~2 times.)
            n = numel(pIdx);
            c = obj.assemblyChunk;
            if n <= c
                A = fun(pIdx,mat);
                return
            end
            perElement = isnumeric(mat) && size(mat,3) == n;
            A = [];
            for s = 1:c:n
                k = s:min(s+c-1,n);
                if perElement
                    Ak = fun(pIdx(k),mat(:,:,k));
                else
                    Ak = fun(pIdx(k),mat);
                end
                if isempty(A), A = Ak; else, A = A + Ak; end
            end
        end

        function attachDOFHandler(obj,DOFHandler)
            obj.DOFs = DOFHandler;
		end
		
		function checkBoundaries(obj,quiet)
            if size(obj.geometry.bd,1) == 3
                nB= size(obj.geometry.bd,2);
                for i=1:nB
                    for j = i+1:nB
						if ~isempty(obj.geometry.bd{3,i}) && ~isempty(obj.geometry.bd{3,j})
							[~,idx] = ismember(obj.geometry.bd{3,i}.nodes,obj.geometry.bd{3,j}.nodes);
							if sum(idx) > 0
								if quiet == 2
									while true
										display('Please choose which to keep!',...
											'2 Boundaries contain the same nodes!');
										x = input(['[1] ',obj.geometry.bd{3,i}.name,' or [2] ',obj.geometry.bd{3,j}.name,'\n>> ']);
										if x == 1
											obj.geometry.bd{3,j}.nodes(nonzeros(idx)) = [];
											[~,idx] = ismember(obj.geometry.bd{3,i}.edges,obj.geometry.bd{3,j}.edges);
											obj.geometry.bd{3,j}.edges(nonzeros(idx)) = [];
											%obj.geometry.bd{3,j}.value(nonzeros(idx)) = [];
											[~,idx] = ismember(obj.geometry.bd{3,i}.faces,obj.geometry.bd{3,j}.faces);
											obj.geometry.bd{3,j}.faces(nonzeros(idx)) = [];
											break;
										elseif x == 2
											obj.geometry.bd{3,i}.nodes(nonzeros(idx)) = [];
											[~,idx] = ismember(obj.geometry.bd{3,i}.edges,obj.geometry.bd{3,j}.edges);
											obj.geometry.bd{3,i}.edges(nonzeros(idx)) = [];
											%obj.geometry.bd{3,j}.value(nonzeros(idx)) = [];
											[~,idx] = ismember(obj.geometry.bd{3,i}.faces,obj.geometry.bd{3,j}.faces);
											obj.geometry.bd{3,i}.faces(nonzeros(idx)) = [];
											break;
										else
											disp('Invalid Input!')
										end
									end
								else
									obj.geometry.bd{3,j}.nodes(nonzeros(idx)) = [];
								end
							end
						end
                    end
                end
            end
		end
        
        function du = gradCell(obj,u)
            % Computes the gradient inside each cell. No patch recovery
			refTet = obj.geometry.refTet;
            DinvT = obj.geometry.DinvT;
            uElem = u(obj.DOFs.el2DOF(:,:));
            uElem = reshape(uElem',size(uElem,2),1,[]);
			if obj.element.dim == 2
				dphi(:,:,1) = obj.element.dPhi{1}(1/3,1/3);
				dphi(:,:,2) = obj.element.dPhi{2}(1/3,1/3);
				du = pagemtimes(pagemtimes(DinvT,dphi(:,:,refTet)),uElem);
			elseif obj.element.dim == 3
				dphi(:,:,1) = obj.element.dPhi{1}(1/4,1/4,1/4);
				dphi(:,:,2) = obj.element.dPhi{2}(1/4,1/4,1/4);
				du = pagemtimes(pagemtimes(DinvT,dphi(:,:,refTet)),uElem);
			end
        end
        
        
        function grad=gradu(obj,u)
        % GRADU computes the gradient at DOFs. 
        %
        % grad=gradu(u) computes the gradient grad of the FEM solution u at
        % DOFs. grad is a Ndofs by Nd matrix, where Ndofs are the number of
        % nodes and Nd the dimension of the spatial space.
        %
        
            switch obj.element.shape_function
                case 'H1'
                        
                    switch obj.geometry.dim
                        case 2
                            %d denotes the dimension of polynomial space,
                            %i.e. for P1 elements it is 1
                            d = 1;             % degree of finite element space
                            m = (d+2)*(d+3)/2; % polynomial degree of approximant
                            n = (d+3)*(d+4)/2; % degree of point considering for least-squares

                            co   = double(squeeze(obj.geometry.co))'; 
                            Ndof = size(co,1);
                            
                            if Ndof < n
                                    error('Mesh too coarse to compute gradient. Please refine it.');
                            end
                            
                            co_idx = knnsearch(co,co,'K',n)';

                            ui1 = reshape(u (co_idx(:)  ),n,1,[]);
                            xi1 = reshape(co(co_idx(:),1),n,1,[]);
                            yi1 = reshape(co(co_idx(:),2),n,1,[]);

                            ui = u(co_idx(:));
                            xi = squeeze(obj.geometry.co(1, :, co_idx(:)));
                            yi = squeeze(obj.geometry.co(2, :, co_idx(:)));

                            clear co_idx;

                            A = [ones(n*Ndof,1), xi, yi, xi.*yi, xi.^2, yi.^2];                    

                            ui_r = cellfun(@(A,b) A\b      , ...
                                            mat2cell(A ,n*ones(1,Ndof),m), ...
                                            mat2cell(ui,n*ones(1,Ndof),1), ...
                                            'UniformOutput',false);
                            ui_r = reshape(cell2mat(ui_r),m,[])';

                            x = co(:,1);
                            y = co(:,2);

                            grad = [ ui_r(:,2)+ui_r(:,4).*y+2*ui_r(:,5).*x, ...
                                     ui_r(:,3)+ui_r(:,4).*x+2*ui_r(:,6).*y ];

                              % DEVELOPMENT statistics and machine learning toolbox
%                             co = double(squeeze(obj.geometry.co))';
%                             I  = obj.geometry.el(:,[1 1 2 2 3 3]);
%                             J  = obj.geometry.el(:,[2 3 1 3 1 2]);
%                             h  = co(I,:)-co(J,:); h=sqrt(dot(h,h,2));
%                             h  = accumarray(I(:),h(:),[],@max,[],true);
%                             ht = full(h);
%                             Di = pdist2(co,co);
%                             while 1
%                                 idx  = Di<repmat(ht,1,Ndof);
%                                 idx2 = sum(idx,2)<m;
%                                 if ~any(idx2)
%                                     break;
%                                 end
%                                 ht(idx2)=ht(idx2)+h(idx2);
%                             end
%                             clear I J idx2 maxR maxRt D;
% 
%                             ui_r=zeros(6,1,Ndof);
%                             for i=1:Ndof
%                                 idxl        = find(idx(i,:));
%                                 ui          = u (idxl  );
%                                 xi          = co(idxl,1);
%                                 yi          = co(idxl,2);
%                                 A           = [ones(numel(idxl),1), xi, yi, xi.*yi, xi.^2, yi.^2];
%                                 ui_r(:,:,i) = A\ui;
%                             end
% 
%                             x = obj.geometry.co(1,:,:);
%                             y = obj.geometry.co(2,:,:);
% 
%                             grad = [ ui_r(2,:,:)+ui_r(4,:,:).*y+2*ui_r(5,:,:).*x, ...
%                                  ui_r(3,:,:)+ui_r(4,:,:).*x+2*ui_r(6,:,:).*y ];
% 
%                             grad = reshape(grad,2,[])';                                 

                            
                    
                    
                        case 3
                            % d denotes the dimension of polynomial space,
                            % i.e. for P1 elements it is 1
                            %d=1;             % degree of finite element space
                            %m=(d+2)*(d+3)/2; % polynomial degree of approximant
                            m = 10;
                            % n=(d+3)*(d+4)/2; % degree of point considering for least-squares
                            n = 3*m;

                            % statistics and machine learning toolbox

                            co_idx = knnsearch(squeeze(obj.geometry.co)',squeeze(obj.geometry.co)','K',n)';

                            %                 ui = u (co_idx(:)  );
                            %                 xi = co(co_idx(:),1);
                            %                 yi = co(co_idx(:),2);
                            %                 zi = co(co_idx(:),3);

                            ui = u(co_idx(:));
                            xi = squeeze(obj.geometry.co(1, :, co_idx(:)));
                            yi = squeeze(obj.geometry.co(2, :, co_idx(:)));
                            zi = squeeze(obj.geometry.co(3, :, co_idx(:)));

                            clear co_idx;

                            Ndof = size(obj.geometry.co,3);

                            A = [ones(Ndof*n,1), xi, yi, zi, xi.*yi ,xi.*zi, yi.*zi, xi.^2, yi.^2, zi.^2];

                            ui_r = cellfun(@(A,b) A\b      , ...
                                mat2cell(A ,n*ones(1,Ndof),m), ...
                                mat2cell(ui,n*ones(1,Ndof),1), ...
                                'UniformOutput',false);
                            ui_r = reshape(cell2mat(ui_r),m,[])';

                            clear A;

                            x = squeeze(obj.geometry.co(1, :,:));
                            y = squeeze(obj.geometry.co(2,:, :));
                            z = squeeze(obj.geometry.co(3, :,:));

                            grad = [ ui_r(:,2)+ui_r(:,5).*y+ui_r(:,6).*z+2*ui_r(:, 8).*x, ...
                                     ui_r(:,3)+ui_r(:,5).*x+ui_r(:,7).*z+2*ui_r(:, 9).*y,...
                                     ui_r(:,4)+ui_r(:,6).*x+ui_r(:,7).*y+2*ui_r(:,10).*z];

                    end
%                     otherwise
%                         error('ofem:elliptic:NotSupported',...
%                               'Reconstruction of P2 gradient not implemented, yet!');
            end
        end
  
        function solve(obj)
            %DOFs = obj.DOFs.getDOFs;
            dofs = obj.DOFs.freeDOFs;
            obj.u(dofs) = obj.S(dofs,dofs)\obj.b(dofs);
            
        end
        
        function parabolic_solve(obj, init, timesteps, deltaTime,folder,file)
            % Implicit Euler for M du/dt + (S + M_robin) u = b with the
            % (time independent) Dirichlet values set by assemble.
            %
            % init.value: a constant initial value (for H1 elements only
            % the vertex DOFs are set, which represents the constant
            % exactly in the hierarchical basis) or a full DOF vector.
            % Exports u after every step to [pwd,folder]/[file,i].
            free = obj.DOFs.freeDOFs;
            dt = deltaTime;

            % Dirichlet values and right hand side without the lifting
            % -(S+M+D+M_robin)*uD that assemble subtracted
            uD = obj.u;
            uD(free) = 0;
            K = obj.S+obj.M+obj.D+obj.M_robin;
            bNat = obj.b + K*uD;

            % initial value
            u = uD;
            if isscalar(init.value)
                if ~isempty(obj.DOFs.n2DOF)
                    vertex = false(size(u));
                    vertex(obj.DOFs.n2DOF(:)) = true;
                    isFree = false(size(u));
                    isFree(free) = true;
                    u(vertex & isFree) = init.value;
                else
                    u(free) = init.value;
                end
            else
                u(free) = init.value(free);
            end
            obj.u = u;

            A = obj.S + obj.M/dt + obj.M_robin;
            Aff = A(free,free);
            try
                dA = decomposition((Aff+Aff')/2,'chol');
            catch
                dA = decomposition(Aff,'lu');
            end
            obj.geometry.export_UCD([pwd, folder],[file,num2str(0)], {'T', full(obj.u), ''});
            for i = 1:timesteps
                rhs = bNat + obj.M/dt*obj.u - A*uD;
                obj.u = uD;
                obj.u(free) = dA\rhs(free);
                obj.geometry.export_UCD([pwd, folder],[file,num2str(i)], {'T', full(obj.u), ''});
            end
        end
        
    end
end


















