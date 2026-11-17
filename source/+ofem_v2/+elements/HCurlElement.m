classdef HCurlElement < ofem_v2.elements.Finite_Elements & handle
    %NEDELECZAGL Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        shape_function;
        gausianWeight;
        gausianLength;
        gradShapeFunctions;
        pointwiseProductShapeFunctions;
        dim;
        degree;
        degreeMass;
        degreeStiff;
        nodeDOFs = 0;
        edgeDOFs = 0;
        faceDOFs = 0;
        interiorDOFs = 0;
        DOFsPerElement;
        Stiffness_matrix;
        Damping_matrix;
        Mass_matrix;
        Load_vector;
        N;
        curlN;
        condense = false;
    end

    methods
        function obj = HCurlElement(dim, deg)
            %NEDELEC Construct an instance of this class
            %   Detailed explanation goes here
            obj.dim = dim;
            obj.degree = deg;
            if deg == 0; deg = 1; end
            % exact for constant materials on affine elements
            obj.degreeMass = deg*2;
            obj.degreeStiff = max(deg-1,0)*2;
        end

        function [N,curlN] = computeBasis(obj)
            switch obj.dim
                case 2
                    syms u v x t c
                    dr = [u,v];

                    l1 = 1-u-v;
                    l2 = u;
                    l3 = v;

                    l = [l1,l2,l3];

                    E{1} = l(2)*gradient(l(1),dr)-l(1)*gradient(l(2),dr);
                    E{2} = l(3)*gradient(l(1),dr)-l(1)*gradient(l(3),dr);
                    E{3} = l(3)*gradient(l(2),dr)-l(2)*gradient(l(3),dr);

                    E = [E{1},E{2},E{3}];

                    leg = ofem_v2.tools.LegPoly;
                    leg.computePolynomials(obj.degree+2);

                    I = [];

                    ke = [1,2;1,3;2,3];

                    for i = 0:obj.degree-1
                        for k = 1:3
                            E = [E,gradient(subs(leg.legIntS(i+2),[x,t],[l(ke(k,1))-l(ke(k,2)),l(ke(k,1))+l(ke(k,2))]),dr)];
                        end
                    end

                    for i = 0:obj.degree-2
                        for j = 0:obj.degree-i-2
                            m = subs(leg.legIntS(i+2),[x,t],[l(2)-l(1),l(1)+l(2)]);
                            n = l(3)*subs(leg.leg(j+1),[x],[2*l(3)-1]);
                            I = [I,gradient(m,dr)*n+m*gradient(n,dr)];
                            I = [I,gradient(m,dr)*n-m*gradient(n,dr)];
                            if i==0
                                I = [I,(gradient(l(1),dr)*l(2)-l(1)*gradient(l(2),dr))*n];
                            end
                        end
                    end

                    E = reshape([E,E],obj.dim,size(E,2),2);
                    I = reshape([I,I],obj.dim,size(I,2),2);
                    N = [E,I];
                    N = simplify(N);
                    obj.DOFsPerElement = size(N,2);
                    %curlN = N;
                    for i = 1:size(N,2)
                        curlN(:,i,1) = curl([N(:,i,1);sym(0)],[dr,c]);
                        curlN(:,i,2) = curl([N(:,i,2);sym(0)],[dr,c]);
                    end

                    curlFunc{1} = matlabFunction(curlN(3,:,1),'vars',[u,v]);
                    curlFunc{2} = matlabFunction(curlN(3,:,2),'vars',[u,v]);
                    NFunc{1} = matlabFunction(N(:,:,1),'vars',[u,v]);
                    NFunc{2} = matlabFunction(N(:,:,2),'vars',[u,v]);

                    obj.edgeDOFs = obj.degree+1;
                    obj.interiorDOFs = size(I,2);
                    N = NFunc;
                    curlN = curlFunc;
                    obj.N = NFunc;
                    obj.curlN = curlFunc;


                case 3
                    %% First compute the zeroth order Nedelec elements to construct everything
                    syms u v w x t c;
                    dr = [u,v,w];

                    %% H1 Basis
                    l1 = 1-u-v-w;
                    l2 = u;
                    l3 = v;
                    l4 = w;

                    l = [l1,l2,l3,l4];


                    %% Zeroth order H curl base
                    % edge 1 2
                    ke = [1,2;1,3;1,4;2,3;2,4;3,4];

                    E{1} = l(2)*gradient(l(1),dr)-l(1)*gradient(l(2),dr);
                    E{2} = l(3)*gradient(l(1),dr)-l(1)*gradient(l(3),dr);
                    E{3} = l(4)*gradient(l(1),dr)-l(1)*gradient(l(4),dr);
                    E{4} = l(3)*gradient(l(2),dr)-l(2)*gradient(l(3),dr);
                    E{5} = l(2)*gradient(l(3),dr)-l(3)*gradient(l(2),dr);
                    E{6} = l(4)*gradient(l(2),dr)-l(2)*gradient(l(4),dr);
                    E{7} = l(4)*gradient(l(3),dr)-l(3)*gradient(l(4),dr);

                    E1 = [E{1},E{2},E{3},E{4},E{6},E{7}];
                    E2 = [E{1},E{2},E{3},E{5},E{6},E{7}];

                    leg = ofem_v2.tools.LegPoly;
                    leg.computePolynomials(obj.degree+2);
                    F1 = {};
                    F2 = {};
                    %                     F{1,4,2} = [];
                    I = {};

                    for i = 0:obj.degree-1
                        for k = 1:6
                            E1 = [E1,gradient(subs(leg.legIntS(i+2),[x,t],[l(ke(k,1))-l(ke(k,2)),l(ke(k,1))+l(ke(k,2))]),dr)];  %Gradient Field
                            if k == 4
                                E2 = [E2,gradient(subs(leg.legIntS(i+2),[x,t],[l(ke(k,2))-l(ke(k,1)),l(ke(k,1))+l(ke(k,2))]),dr)];
                            else
                                E2 = [E2,gradient(subs(leg.legIntS(i+2),[x,t],[l(ke(k,1))-l(ke(k,2)),l(ke(k,1))+l(ke(k,2))]),dr)];
                            end
                        end
                    end

                    kf1 = [1,2,3;1,2,4;1,3,4;2,3,4];
                    kf2 = [1,3,2;1,2,4;1,3,4;3,2,4];

                    for k = 1:4
                        for i =0:obj.degree-2
                            for j = 0:obj.degree-i-2
                                m1 = subs(leg.legIntS(i+2),[x,t],[l(kf1(k,1))-l(kf1(k,2)),l(kf1(k,1))+l(kf1(k,2))]);
                                n1 = l(kf1(k,3))*subs(leg.legS(j+1),[x,t],[l(kf1(k,3))-l(kf1(k,1))-l(kf1(k,2)),l(kf1(k,3))+l(kf1(k,1))+l(kf1(k,2))]);
                                m2 = subs(leg.legIntS(i+2),[x,t],[l(kf2(k,1))-l(kf2(k,2)),l(kf2(k,1))+l(kf2(k,2))]);
                                n2 = l(kf2(k,3))*subs(leg.legS(j+1),[x,t],[l(kf2(k,3))-l(kf2(k,1))-l(kf2(k,2)),l(kf2(k,3))+l(kf2(k,1))+l(kf2(k,2))]);
                                F1 = [F1,gradient(m1*n1,dr)];     % Gradient Field
                                F1 = [F1,gradient(n1,dr)*m1-n1*gradient(m1,dr)];
                                F2 = [F2,gradient(m2*n2,dr)];     % Gradient Field
                                F2 = [F2,gradient(n2,dr)*m2-n2*gradient(m2,dr)];
                                if i==0
                                    F1 = [F1,(gradient(l(kf1(k,1)),dr)*l(kf1(k,2))-l(kf1(k,1))*gradient(l(kf1(k,2)),dr))*n1];
                                    F2 = [F2,(gradient(l(kf2(k,1)),dr)*l(kf2(k,2))-l(kf2(k,1))*gradient(l(kf2(k,2)),dr))*n2];
                                end
                            end
                        end
                    end
                    for i = 0:obj.degree-3
                        for j = 0:obj.degree-3-i
                            for k = 0:obj.degree-3-i-j
                                m = subs(leg.legIntS(i+2),[x,t],[l(1)-l(2),l(1)+l(2)]);
                                n = l(3)*subs(leg.legS(j+1),[x,t],[2*l(3)-(1-l(4)),1-l(4)]);
                                p = l(4)*subs(leg.leg(k+1),x,2*l(4)-1);
                                I = [I,gradient(m*n*p,dr)];     %% Gradient Field
                                I = [I,gradient(m,dr)*n*p-m*gradient(n,dr)*p+m*n*gradient(p,dr)];
                                I = [I,gradient(m,dr)*n*p+m*gradient(n,dr)*p-m*n*gradient(p,dr)];
                                if i==0
                                    I = [I,(l(2)*gradient(l(1),dr)-l(1)*gradient(l(2),dr))*n*p];
                                end
                            end
                        end
                    end
                    I2 = I;
                    E = reshape([E1,E2],obj.dim,size(E1,2),2);
                    F = reshape([F1,F2],obj.dim,size(F1,2),2);
                    I = reshape([I,I2],obj.dim,size(I,2),2);
                    N = [E,F,I];
                    N = simplify(N);
                    obj.DOFsPerElement = size(N,2);
                    curlN = N;
                    for i = 1:size(curlN,2)
                        curlN(:,i,1) = curl(curlN(:,i,1),dr);
                        curlN(:,i,2) = curl(curlN(:,i,2),dr);
                    end

                    curlFunc{1} = matlabFunction(curlN(:,:,1),'vars',[u,v,w]);
                    curlFunc{2} = matlabFunction(curlN(:,:,2),'vars',[u,v,w]);
                    NFunc{1} = matlabFunction(N(:,:,1),'vars',[u,v,w]);
                    NFunc{2} = matlabFunction(N(:,:,2),'vars',[u,v,w]);

                    obj.faceDOFs = size(F,2)/4;
                    obj.edgeDOFs = obj.degree+1;
                    obj.interiorDOFs = size(I,2);
                    N = NFunc;
                    curlN = curlFunc;
                    obj.N = NFunc;
                    obj.curlN = curlFunc;
            end
        end



        function S = assembleStiffness(obj,phys,pIdx,mat)
            % int curl(u)' * mat * curl(v)
            refTet = phys.geometry.refTet(pIdx);
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeStiff);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            Dk = phys.geometry.Dk(:,:,pIdx);
            if obj.dim == 2
                Dk = [];        % the 2D curl is a scalar: curl u = curl_ref/detD
            end

            if isa(mat,'function_handle')
                Ne = length(pIdx);
                S = zeros(obj.DOFsPerElement,obj.DOFsPerElement,Ne);
                X = quadPoints(phys,pIdx,l);
                for q = 1:size(w,1)
                    dphi = evalRef(obj.curlN,l(:,q),refTet);
                    if ~isempty(Dk)
                        dphi = pagemtimes(Dk,dphi);
                    end
                    S = S+w(q)*pagemtimes(dphi,'transpose',materialTimes(mat(X{q}),dphi),'none');
                end
                S = S./abs(detD);
            else
                % constant material: reference tensors
                S = ofem_v2.tools.refTensorMatrices(obj.curlN,w,l,refTet, ...
                    ofem_v2.tools.geometryFactor(Dk,mat,1./abs(detD)));
            end
            S = toSparse(phys,dofs,S);
        end

        function M = assembleMass(obj,phys,pIdx,mat)
            % int u' * mat * v
            refTet = phys.geometry.refTet(pIdx);
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            DinvT = phys.geometry.DinvT(:,:,pIdx);

            if isa(mat,'function_handle')
                Ne = length(pIdx);
                M = zeros(obj.DOFsPerElement,obj.DOFsPerElement,Ne);
                X = quadPoints(phys,pIdx,l);
                for q = 1:size(w,1)
                    phi = pagemtimes(DinvT,evalRef(obj.N,l(:,q),refTet));
                    M = M+w(q)*pagemtimes(phi,'transpose',materialTimes(mat(X{q}),phi),'none');
                end
                M = M.*abs(detD);
            else
                % constant material: reference tensors
                M = ofem_v2.tools.refTensorMatrices(obj.N,w,l,refTet, ...
                    ofem_v2.tools.geometryFactor(DinvT,mat,abs(detD)));
            end
            M = toSparse(phys,dofs,M);
        end

        function D = assembleDamping(obj,phys,pIdx,mat) %#ok<INUSD>
            D = [];
        end

        function D = assembleVelocityTerm(obj,phys,pIdx,mat,v)
            % int u' * mat * (v x curl(w)), v a constant 3 x 1 vector, mat
            % scalar, 3 x 3 or a function handle of the coordinates. In 2D
            % u is extended by a zero third component and curl(w) = (0,0,curl).
            refTet = phys.geometry.refTet(pIdx);
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            DinvT = phys.geometry.DinvT(:,:,pIdx);
            Dk = phys.geometry.Dk(:,:,pIdx);
            Ns = obj.DOFsPerElement;
            Ne = length(pIdx);
            if isa(mat,'function_handle')
                X = quadPoints(phys,pIdx,l);
            end
            v = repmat(v(:),1,Ns,Ne);

            D = zeros(Ns,Ns,Ne);
            for q = 1:size(w,1)
                phi = pagemtimes(DinvT,evalRef(obj.N,l(:,q),refTet));
                c = evalRef(obj.curlN,l(:,q),refTet);
                if obj.dim == 3
                    c = pagemtimes(Dk,c)./detD;
                else
                    phi(3,:,:) = 0;
                    c = [zeros(2,Ns,Ne); c./detD];
                end
                vc = cross(v,c,1);
                if isa(mat,'function_handle')
                    vc = materialTimes(mat(X{q}),vc);
                else
                    vc = materialTimes(mat,vc);
                end
                D = D+w(q)*pagemtimes(phi,'transpose',vc,'none');
            end
            D = toSparse(phys,dofs,D.*abs(detD));
        end

        function F = volumeForce(obj,phys,f,pIdx)
            % int f . v, f a constant vector, a function handle of the
            % coordinates, values per element (dim x 1 x Ne) or at the
            % quadrature points (dim x Nq x Ne)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            DinvT = phys.geometry.DinvT(:,:,pIdx);
            refTet = phys.geometry.refTet(pIdx);
            Ne = length(pIdx);
            F = zeros(1,obj.DOFsPerElement,Ne);
            X = {};
            if isa(f,'function_handle')
                X = quadPoints(phys,pIdx,l);
            elseif size(f,3) > 1
                f = double(f(:,:,pIdx));
            end
            for q = 1:size(w,1)
                phi = pagemtimes(DinvT,evalRef(obj.N,l(:,q),refTet));
                F = F+w(q)*pagemtimes(loadValues(f,X,q),'transpose',phi,'none');
            end
            F = F.*abs(detD);
            F = sparse(reshape(dofs',[],1),1,F(:),phys.DOFs.Nd,1);
        end

        function F = derivativeForce(obj,phys,f,pIdx)
            % int f . curl(v), f as in volumeForce (a scalar field in 2D)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            Dk = phys.geometry.Dk(:,:,pIdx);
            refTet = phys.geometry.refTet(pIdx);
            Ne = length(pIdx);
            F = zeros(1,obj.DOFsPerElement,Ne);
            X = {};
            if isa(f,'function_handle')
                X = quadPoints(phys,pIdx,l);
            elseif size(f,3) > 1
                f = double(f(:,:,pIdx));
            end
            for q = 1:size(w,1)
                c = evalRef(obj.curlN,l(:,q),refTet);
                if obj.dim == 3
                    c = pagemtimes(Dk,c);
                end
                F = F+w(q)*pagemtimes(loadValues(f,X,q),'transpose',c,'none');
            end
            % curl v = (Dk) curl_ref/detD, dx = |detD| dx_ref
            F = F.*sign(detD);
            F = sparse(reshape(dofs',[],1),1,F(:),phys.DOFs.Nd,1);
        end

        function U = el2edges(obj,phys,u)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF;
            detD = phys.geometry.detD;
            DinvT = phys.geometry.DinvT;

            Ns = obj.DOFsPerElement;
            Nq = size(w   ,1);
            Ne = phys.geometry.Nint;
            Nl = size(phys.geometry.el,2);

            U = ofem_v2.tools.matrixarray(zeros(1,Ns,Ne));

            refTet = phys.geometry.refTet;

            for q=1:Nq
                phi(:,:,1) = obj.N{1}(l(1,q),l(2,q),l(3,q));
                phi(:,:,2) = obj.N{2}(l(1,q),l(2,q),l(3,q));
                phi =  DinvT*ofem_v2.tools.matrixarray(phi(:,:,refTet));
                U = U + w(q)*(u'*phi);
            end

            U = U'*abs(detD);
            I = dofs';

            U = sparse(I(:),1,U(:),phys.DOFs.Nd,1);
            U = full(U);
        end

        function [X,U] = curl(obj,phys,u)
            % First compute the global evaluation points
            [~,l] = ofem_v2.tools.gaussSimplex(phys.element.dim,1);%phys.element.degreeMass);
            Dk = phys.geometry.Dk;
            detD = phys.geometry.detD;
            Ne = phys.geometry.Nint;
            Nl = 4;
            Np = size(l,2);
            elco = reshape(phys.geometry.co(:,:,phys.geometry.el(:,1:Nl)'),[],Nl,Ne);
            X = zeros(3,Ne,Np);
            for i = 1:Np
                % local vertex 1 has barycentric 1-sum(l)
                X(:,:,i) = reshape(pagemtimes(elco,[1-sum(l(:,i));l(:,i)]),3,[]);
            end
            elu = ofem_v2.tools.matrixarray(reshape(u(phys.DOFs.el2DOF(:,:)'),[],1,Ne));
            for q = 1:Np
                dphi(:,:,1) = phys.element.curlN{1}(l(1,q),l(2,q),l(3,q));
                dphi(:,:,2) = phys.element.curlN{2}(l(1,q),l(2,q),l(3,q));
                dphi = ofem_v2.tools.matrixarray(dphi(:,:,phys.geometry.refTet));
                U(:,:,q) = reshape((Dk*dphi)*elu*ofem_v2.tools.matrixarray(1/(detD)),3,[]);
            end
        end

        function Y_g = AMS(obj,phys)
            Ned = phys.geometry.Ned;
            Nco = phys.geometry.Nco;
            Nfa = phys.geometry.Nfa;
            Nint = phys.geometry.Nint;

            Y_N0 = ones(2*Ned,1);
            Y_N0(1:Ned) = -1;

            I = repmat(1:Ned,1,2);
            J = phys.geometry.ed(:);

            Y_N0 = sparse(I(:),J,Y_N0,Ned,Nco);

            Y_E = speye(Ned*obj.degree);

            if obj.dim == 3 && obj.faceDOFs>0
                sten1 = sparse(diag([1,0,0]));
                sten2 = sparse(diag([1,0]));

                eye1 = speye(obj.degree-1);
                eye2 = speye(sum(1:obj.degree-2));

                stencil = blkdiag(kron(eye1,sten1),kron(eye2,sten2));

                Y_F = kron(speye(Nfa),stencil);
            else
                Y_F = [];
            end

            if obj.dim == 3 && obj.degree > 2 && obj.interiorDOFs > 0
                sten1 = sparse(diag([1,0,0,0]));
                sten2 = sparse(diag([1,0,0]));

                eye1 = speye(sum(1:obj.degree-2));
                eye2 = speye(sum(sum(triu(repelem(1:obj.degree-3,obj.degree-3,1)'))));

                stencil = blkdiag(kron(eye1,sten1),kron(eye2,sten2));

                Y_I = kron(speye(Nint),stencil);

            elseif obj.dim == 2 && obj.interiorDOFs > 0
                sten1 = sparse(diag([1,0,0]));
                sten2 = sparse(diag([1,0]));

                eye1 = speye(sum(1:obj.degree-1));
                eye2 = speye(sum(sum(triu(repelem(1:obj.degree-2,obj.degree-2,1)'))));

                stencil = blkdiag(kron(eye1,sten1),kron(eye2,sten2));

                Y_I = kron(speye(Nint),stencil);
            else
                Y_I = [];
            end

            Y_g = blkdiag(Y_N0,Y_E,Y_F,Y_I);
        end

    end

end

function B = materialTimes(A,P)
    % material values A (1 x 1 x Ne scalar or d x d x Ne tensor) times the
    % basis values P (d x Ns x Ne)
    A = double(A);
    if size(A,1) == 1 && size(A,2) == 1
        B = A.*P;
    else
        B = pagemtimes(A,P);
    end
end

function P = evalRef(fun,lq,refTet)
    % basis (or curl) at the reference point lq for both reference
    % elements, selected per element
    c = num2cell(lq);
    P = cat(3,fun{1}(c{:}),fun{2}(c{:}));
    P = P(:,:,refTet);
end

function X = quadPoints(phys,pIdx,l)
    % physical coordinates of the quadrature points, X{q}: dim x 1 x Ne
    g = phys.geometry;
    Nl = size(l,1)+1;
    elco = reshape(g.co(:,:,g.el(pIdx,1:Nl)'),[],Nl,numel(pIdx));
    X = cell(1,size(l,2));
    for q = 1:size(l,2)
        % local vertex 1 has barycentric 1-sum(l)
        X{q} = pagemtimes(elco,[1-sum(l(:,q));l(:,q)]);
    end
end

function A = loadValues(f,X,q)
    % load values at quadrature point q: c x 1 x Ne or c x 1
    if isa(f,'function_handle')
        A = double(f(X{q}));
    elseif size(f,3) > 1 && size(f,2) > 1
        A = f(:,q,:);
    else
        A = double(f);
    end
end

function A = toSparse(phys,dofs,A)
    Ns = size(A,1);
    I = repmat(dofs,1,Ns)';
    J = repelem(dofs',Ns,1);
    A = sparse(I(:),J(:),A(:),phys.DOFs.Nd,phys.DOFs.Nd);
end
