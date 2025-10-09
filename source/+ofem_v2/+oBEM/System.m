classdef System < handle
    %SYSTEM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mesh;
        fe;
        assemble_sp;
        assemble_dp;
        normals;
        distances;
        SP;
        DP;
        M;
        gamma;
        index;
        name;
        Dk, DinvT, detD;
        r;
    end

    methods
        function obj = System(mesh, fe, ID, has_sp, has_dp)
            obj.mesh = mesh;
            obj.fe = fe;
            obj.assemble_sp = has_sp;
            obj.assemble_dp = has_dp;
            n = size(mesh.bd, 2);
            for i = 1:n
                if isequal(lower(mesh.bd{1,i}),lower(ID))
                    obj.gamma = mesh.bd{2,i};
                    obj.index = i;
                    obj.name = ID;
                    break;
                end
            end

        end

        function assemble(obj)
            obj.computeGeometricData;
            if obj.assemble_sp
                obj.SP = obj.assembleSingle;
            end
            if obj.assemble_dp
                obj.DP = obj.assembleDouble;
            end
            obj.M = obj.lumpedMass;

        end

        function computeGeometricData(obj)
            e12 = obj.mesh.co(:,:,obj.gamma(:,2))-obj.mesh.co(:,:,obj.gamma(:,1));
            e13 = obj.mesh.co(:,:,obj.gamma(:,3))-obj.mesh.co(:,:,obj.gamma(:,1));
            e12 = double(e12);
            e13 = double(e13);

            obj.Dk = [e12,e13];
            obj.DinvT = pagemtimes(obj.Dk,'transpose',obj.Dk,'none');
            obj.DinvT = pageinv(obj.DinvT);
            obj.normals = cross(e12,e13,1);
            obj.detD = vecnorm(obj.normals);
            obj.normals = pagemtimes(obj.normals,pageinv(obj.detD));

            obj.r = double(obj.mesh.co(:,:,obj.gamma'));
            obj.r = reshape(obj.r,3,3,[]);
        end

        function S = assembleSingle(obj)

            [S,I,J] = obj.single;

            S = sparse(I,J,S,obj.mesh.Nco,obj.mesh.Nco);
            S = S+S';

            [Ss,Is,Js] = obj.singleSingular2;

            Ss = sparse(Js,Is,Ss,obj.mesh.Nco,obj.mesh.Nco);

            S = S + Ss;
        end

        function [S,I,J] = single(obj)
            [wi,li] = ofem_v2.tools.gaussSimplex(2,2);
            [wj,lj] = ofem_v2.tools.gaussSimplex(2,3);

            li(3,:) = 1-sum(li,1);
            lj(3,:) = 1-sum(lj,1);

            % Create reduced system -> 1,2 1,3 2,3 1,4 2,4 3,4...
            % This way we only compute the upper triangular matrix
            % Singular diagonal is handled seperately
            idx = ones(size(obj.gamma,1),size(obj.gamma,1)).*[1:size(obj.gamma,1)];
            idxi = tril(idx(1:end-1,1:end-1))';
            idxi(idxi==0) = [];
            vi = obj.r(:,:,idxi);
            idxj = triu(idx,1);
            idxj(idxj==0) = [];
            vj = obj.r(:,:,idxj);

            S = zeros(3,3,size(idxi,1));

            for i=1:size(wi,1)
                ri = pagemtimes(vi,li([3,1,2],i));
                cnt = ones(size(li(1:2,i),1),1);
				lTemp = mat2cell(li(1:2,i),cnt);
                phii(:,:) = obj.fe.phi{1}(lTemp{:});
                for j = 1:size(wj,1)
                    rj = pagemtimes(vj,lj([3,1,2],j));
                    cnt = ones(size(lj(1:2,j),1),1);
					lTemp = mat2cell(lj(1:2,j),cnt);
                    phij(:,:) = obj.fe.phi{1}(lTemp{:});

                    rn = pagenorm(ri-rj);
                    G = 1/(4*pi)*pageinv(rn);

                    w_tot = wi(i)*wj(j)*pagemtimes(obj.detD(:,:,idxi),obj.detD(:,:,idxj));
                    cont = pagemtimes(pagemtimes(phii,'transpose',phij,'none'),G);

                    S = S + pagemtimes(w_tot,cont);
                end
            end
            %S = pagemtimes(S,pagemtimes(obj.detD(:,:,idxi),obj.detD(:,:,idxj)));
            %S = pagemtimes(S,'transpose',[1/3,1/3,1/3],'none');

            I = obj.gamma(idxi,:)';
            I = repmat(I,3,1);
            I = I(:);

            J = obj.gamma(idxj,:)';
            J = repelem(J,3,1);
            J = J(:);

            S = S(:);
        end

        function [S,I,J] = singleSingular2(obj)
            [wi,li] = gquts7;%ofem_v2.tools.gaussSimplex(2,6);
            [wj,lj] = gqutm9;%ofem_v2.tools.gaussSimplex(2,6);

            li(3,:) = 1-sum(li,1);
            lj(3,:) = 1-sum(lj,1);

            % Create reduced system -> 1,2 1,3 2,3 1,4 2,4 3,4...
            % This way we only compute the upper triangular matrix
            % Singular diagonal is handled seperately
            vi = obj.r;
            vj = obj.r;

            S = zeros(3,3,size(obj.gamma,1));

            for i=1:size(wi,1)
                ri = pagemtimes(vi,li([3,1,2],i));
                cnt = ones(size(li(1:2,i),1),1);
				lTemp = mat2cell(li(1:2,i),cnt);
                phii(:,:) = obj.fe.phi{1}(lTemp{:});
                for j = 1:size(wj,1)
                    rj = pagemtimes(vj,lj([3,1,2],j));
                    cnt = ones(size(lj(1:2,j),1),1);
					lTemp = mat2cell(lj(1:2,j),cnt);
                    phij(:,:) = obj.fe.phi{1}(lTemp{:});

                    rn = pagenorm(ri-rj);
                    G = 1/(4*pi)*pageinv(rn);

                    w_tot = wi(i)*wj(j);
                    cont = pagemtimes(pagemtimes(phii,'transpose',phij,'none'),G);

                    S = S + pagemtimes(w_tot,cont);
                end
            end
            %S = pagemtimes(S,pagemtimes(obj.detD(:,:,idxi),obj.detD(:,:,idxj)));
            %S = pagemtimes(S,'transpose',[1/3,1/3,1/3],'none');

            S = pagemtimes(S,obj.detD.^2)/4;
            I = repmat(obj.gamma',3,1);
            I = I(:);
            J = repelem(obj.gamma',3,1);
            J = J(:);
            S = S(:);
        end

        function [S,I,J] = singleSingular(obj)
            [wi,li] = ofem_v2.tools.gaussSimplex(1,4);

            S = zeros(3,3,size(obj.r,3));
            AB = pagenorm(obj.Dk(:,2,:)-obj.Dk(:,1,:));
            test = 1;
            for i = 1:size(wi,1)
                for j = 1:size(wi,1)
                    for k = 1:size(wi,1)
                        x = li(i);
                        y = li(j);
                        z = li(k);
                        wiwjwk = wi(i)*wi(j)*wi(k);

                        li = [x*(1-y);x*y];
                        lj = [x*(1-z);x*z];
    
                        cnt = ones(size(li(:),1),1);
				        lTemp = mat2cell(li(:),cnt);
                        phii(:,:) = obj.fe.phi{1}(lTemp{:});
    
                        cnt = ones(size(lj(:),1),1);
				        lTemp = mat2cell(lj(:),cnt);
                        phij(:,:) = obj.fe.phi{1}(lTemp{:});

                        %ri = pagemtimes(obj.r,[1-sum(li,1);li]);
                        %rj = pagemtimes(obj.r,[1-sum(lj,1);lj]);
                        ri = obj.r(:,1,:) + obj.Dk(:,1,:)*li(1) + obj.Dk(:,2,:)*li(2);
                        rj = obj.r(:,1,:) + obj.Dk(:,1,:)*lj(1) + obj.Dk(:,2,:)*lj(2);

                        drtest = pagenorm(rj-ri);

                        % lambda = sqrt(pagemtimes(A,(1-y-z)^2)+ ...
                        %     pagemtimes(B,(y-z)^2)+...
                        %     2*pagemtimes(C,(1-y-z)*(y-z)));
                        lambda = x*abs(y-z)*AB;
    
                        dr = pagemtimes(x,lambda);
                        % dr = lambda;
                         epsdist = sqrt(eps);
                        small = dr < epsdist;
                        G = zeros(1,1,size(obj.gamma,1));
                        % r = pagenorm(ri-rj);
                        G(:,:,~small) = 1/(4*pi)*pageinv(dr);
                        G(:,:,small) = 0;
                        cont = pagemtimes(phii,'transpose',phij,'none');
                        cont = pagemtimes(cont,G);
                        Stest(:,:,test) = wiwjwk*x^2*cont;
                        test = test+1;
                        S = S + wiwjwk*x^2*cont;
                    end
                end
            end
            S = pagemtimes(S,obj.detD.^2);
            I = repmat(obj.gamma',3,1);
            I = I(:);
            J = repelem(obj.gamma',3,1);
            J = J(:);
            S = S(:);

        end

        function D = assembleDouble(obj)
            phi = [1/2,1/2,0; 1/2,0,1/2; 0,1/2,1/2];
            w = [1/3,1/3,1/3];
            factor = 1/(4*pi);

            D = zeros(3,3,size(obj.gamma,1),size(obj.gamma,1));

            for q=1:size(w,2)
                rq = pagemtimes(obj.r,phi(:,q));
                rq = reshape(rq,3,1,1,[]);

                vec = obj.r-rq;
                dist = vecnorm(vec);

                dGdn = -(pagemtimes(vec,"transpose",obj.normals,'none'));
                dGdn = dGdn./(factor*pagetranspose(dist).^3);

                D = D + w(q)*pagemtimes(dGdn,'none',phi(:,q),'transpose');
            end
            D = pagemtimes(D,abs(obj.detD)/2);

            I = repelem(obj.gamma',3,1);
            I = repmat(I,size(obj.gamma,1),1);
            I = I(:);

            J = repmat(obj.gamma',3,1);
            J = J(:);
            J = repmat(J,1,size(obj.gamma,1));
            J = J(:);

            D = sparse(I,J,D(:));
        end

        function M = lumpedMass(obj)
            areas = obj.detD/2;
            areas = pagemtimes(areas,[1/3,1/3,1/3]);
            I = obj.gamma';
            M = sparse(I(:),1,areas(:));
            M = diag(M);
        end
    end
end

