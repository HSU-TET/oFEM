classdef AdditiveJacobi <  handle
    %ADDITIVEJACOBI undefined
    %   undefined

    properties
        prec;                   % Preconditioner matrix
        G;                      % Discrete Gradient
        L;                      % Nodal Matrix
        nodalCorrect = 0;       % 0 or 1 if nodal solve should be used
    end

    methods
        function obj = AdditiveJacobi(dofhandler,A)
            i = size(dofhandler.el2DOF,2);
            j = i;
            k = size(dofhandler.el2DOF,1);
            BJ = zeros(i,j,k);

            for i = 1:k
                BJ(:,:,i) = full(A(dofhandler.el2DOF(i,:),dofhandler.el2DOF(i,:)));
            end

            I = repmat(dofhandler.el2DOF(:,:),1,size(BJ,1))';
            I = I(:);
            J = repelem(dofhandler.el2DOF(:,:)',size(BJ,1),1);
            J = J(:);

            BJ = pageinv(BJ);

            BJ = sparse(I,J,BJ(:),dofhandler.Nd,dofhandler.Nd);
            obj.prec = BJ(dofhandler.freeDOFs,dofhandler.freeDOFs);
        end

        function nodalSolve(obj,G,L,freeDOFs,amsdofs)
            %METHOD1 undefined
            %   undefined
            obj.nodalCorrect = 1;
            obj.G = G(freeDOFs,amsdofs);
            obj.L = decomposition(L(amsdofs,amsdofs),'lu');
        end

        function y = preconditioner(obj,r)
            y_edge = obj.prec*r;

            %r = r-A*y_edge;

            r_nodal = obj.G'*r;
            e_nodal = obj.L\r_nodal;

            y = y_edge + obj.G*e_nodal;
        end
    end
end