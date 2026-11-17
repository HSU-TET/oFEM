classdef HDivElement < ofem_v2.elements.Finite_Elements & handle
    %HDIVELEMENT H(div)-conforming high-order element following Zaglmayr
    %
    %   fe = HDivElement(dim, p) on triangles (dim = 2) or tetrahedra
    %   (dim = 3).
    %
    %   p = 0: lowest-order Raviart-Thomas element (one DOF per facet).
    %   p > 0: full polynomial space P_p^dim (Brezzi-Douglas-Marini).
    %
    %   The basis follows S. Zaglmayr, High Order Finite Element Methods
    %   for Electromagnetic Field Computation, PhD thesis, JKU Linz, 2006:
    %
    %   3D, face f = [f1,f2,f3] in global vertex order,
    %       u_i = L^S_{i+2}(l_f1-l_f2, l_f1+l_f2),
    %       v_j = l_f3 P^S_j(l_f3-l_f1-l_f2, l_f1+l_f2+l_f3)
    %     - Raviart-Thomas: l_f1 Dl_f2 x Dl_f3 + cyclic
    %     - Dv_j x Du_i,               i+j <= p-1   (curl of H(curl) type 2)
    %     - curl(N0_f1f2 v_j),         j   <= p-1   (curl of H(curl) type 3)
    %   3D, cell with u_i = L^S_{i+2}(l1-l2, l1+l2),
    %       v_j = l3 P^S_j(2 l3-(1-l4), 1-l4), w_k = l4 P_k(2 l4-1)
    %     divergence free (curls of H(curl) cell functions)
    %     - Dv_j x D(u_i w_k), Dw_k x D(u_i v_j),   i+j+k <= p-2
    %     - curl(N0_12 v_j w_k),                    j+k   <= p-2
    %     not divergence free
    %     - w_k Du_i x Dv_j,                        i+j+k <= p-2
    %     - l2 w_k Dl1 x Dv_j,                      j+k   <= p-2
    %     - l2 v_j Dl1 x Dl4,                       j     <= p-2
    %   2D: the H(curl) basis of HCurlElement rotated by -90 degrees, so
    %   normal traces are the H(curl) tangential traces and the
    %   divergence is the H(curl) curl.
    %
    %   The functions are mapped with the contravariant Piola transform
    %   u = Dk*u_ref/detD, div u = div u_ref/detD (signed detD).

    properties (Access = public)
        shape_function = 'HDiv';
        dim;
        degree;
        degreeMass;
        degreeStiff;
        nodeDOFs = 0;
        edgeDOFs = 0;
        faceDOFs = 0;
        interiorDOFs = 0;
        DOFsPerElement;
        N;          % basis on the two reference elements, dim x Ns
        divN;       % divergence of the basis, 1 x Ns
    end

    methods
        function obj = HDivElement(dim, deg)
            if ~any(dim == [2 3])
                error('ofem_v2:elements:HDivElement', ...
                    'HDivElement is implemented for dim = 2 and 3.');
            end
            obj.dim = dim;
            obj.degree = deg;
            obj.degreeMass = 2*max(deg,1);
            obj.degreeStiff = 2*max(deg-1,0);
        end

        function [N,divN] = computeBasis(obj)
            switch obj.dim
                case 2
                    [N,divN] = obj.basis2D();
                case 3
                    [N,divN] = obj.basis3D();
            end
            obj.N = N;
            obj.divN = divN;
        end

        function S = assembleStiffness(obj,phys,pIdx,mat)
            % int mat * div(u) * div(v)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeStiff);
            [refTet,dofs,detD,~,X] = obj.elementData(phys,pIdx,l);
            Ns = obj.DOFsPerElement;
            S = zeros(Ns,Ns,numel(pIdx));
            for q = 1:size(w,1)
                dphi = obj.evalRef(obj.divN,l(:,q),refTet);
                S = S + w(q)*pagemtimes(dphi,'transpose',materialTimes(mat,X,q,dphi),'none');
            end
            S = S./abs(detD);
            S = obj.toSparse(phys,dofs,S);
        end

        function M = assembleMass(obj,phys,pIdx,mat)
            % int u' * mat * v
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            [refTet,dofs,detD,Dk,X] = obj.elementData(phys,pIdx,l);
            Ns = obj.DOFsPerElement;
            M = zeros(Ns,Ns,numel(pIdx));
            for q = 1:size(w,1)
                phi = pagemtimes(Dk,obj.evalRef(obj.N,l(:,q),refTet));
                M = M + w(q)*pagemtimes(phi,'transpose',materialTimes(mat,X,q,phi),'none');
            end
            M = M./abs(detD);
            M = obj.toSparse(phys,dofs,M);
        end

        function D = assembleDamping(obj,phys,pIdx,mat) %#ok<INUSD>
            D = [];
        end

        function F = volumeForce(obj,phys,f,pIdx)
            % int f . v, f a dim x 1 vector, a function handle of the
            % coordinates (dim x 1 x Ne) or values per element
            % (dim x 1 x Ne or dim x Nq x Ne at the quadrature points)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            [refTet,dofs,detD,Dk,X] = obj.elementData(phys,pIdx,l);
            Ne = numel(pIdx);
            F = zeros(1,obj.DOFsPerElement,Ne);
            for q = 1:size(w,1)
                phi = pagemtimes(Dk,obj.evalRef(obj.N,l(:,q),refTet));
                if isa(f,'function_handle')
                    fq = double(f(X{q}));
                elseif size(f,3) > 1 && size(f,2) > 1
                    fq = f(:,q,pIdx);
                elseif size(f,3) > 1
                    fq = f(:,:,pIdx);
                else
                    fq = repmat(f(:),1,1,Ne);
                end
                F = F + w(q)*pagemtimes(fq,'transpose',phi,'none');
            end
            F = F.*sign(detD);
            F = sparse(reshape(dofs',[],1),1,F(:),phys.DOFs.Nd,1);
        end

        function [X,U,divU] = evaluate(obj,phys,u,l)
            % Evaluate the discrete field u (and its divergence) at the
            % barycentric points l (dim x Np) of all elements.
            % X, U: dim x Ne x Np, divU: 1 x Ne x Np
            Ne = phys.geometry.Nint;
            pIdx = 1:Ne;
            [refTet,dofs,detD,Dk,Xc] = obj.elementData(phys,pIdx,l);
            elu = reshape(full(u(dofs')),[],1,Ne);
            Np = size(l,2);
            X = zeros(obj.dim,Ne,Np); U = X; divU = zeros(1,Ne,Np);
            for q = 1:Np
                X(:,:,q) = reshape(Xc{q},obj.dim,[]);
                phi = pagemtimes(Dk,obj.evalRef(obj.N,l(:,q),refTet));
                U(:,:,q) = reshape(pagemtimes(phi,elu)./detD,obj.dim,[]);
                dphi = obj.evalRef(obj.divN,l(:,q),refTet);
                divU(:,:,q) = reshape(pagemtimes(dphi,elu)./detD,1,[]);
            end
        end
    end

    methods (Access = private)
        function [N,divN] = basis2D(obj)
            % rotate the H(curl) basis: (n1,n2) -> (n2,-n1)
            hc = ofem_v2.elements.HCurlElement(2,obj.degree);
            [Nc,curlNc] = hc.computeBasis;
            R = [0 1;-1 0];
            N = {@(u,v) R*Nc{1}(u,v), @(u,v) R*Nc{2}(u,v)};
            divN = curlNc;
            obj.edgeDOFs = hc.edgeDOFs;
            obj.interiorDOFs = hc.interiorDOFs;
            obj.DOFsPerElement = hc.DOFsPerElement;
        end

        function [N,divN] = basis3D(obj)
            syms u v w x t
            dr = [u,v,w];
            p = obj.degree;
            l = [1-u-v-w,u,v,w];
            g = cell(1,4);
            for i = 1:4
                g{i} = gradient(l(i),dr);
            end
            grad = @(f) gradient(f,dr);
            rotN = @(F) curl(F,dr);

            leg = ofem_v2.tools.LegPoly;
            leg.computePolynomials(p+2);
            uF = @(i,a,b) subs(leg.legIntS(i+2),[x,t],[l(a)-l(b),l(a)+l(b)]);
            vF = @(j,a,b,c) l(c)*subs(leg.legS(j+1),[x,t],[l(c)-l(a)-l(b),l(a)+l(b)+l(c)]);

            % faces in global vertex order for both reference tetrahedra
            kf = {[1,2,3;1,2,4;1,3,4;2,3,4], [1,3,2;1,2,4;1,3,4;3,2,4]};
            F = cell(1,2);
            for var = 1:2
                F{var} = sym(zeros(3,0));
                for k = 1:4
                    a = kf{var}(k,1); b = kf{var}(k,2); c = kf{var}(k,3);
                    % Raviart-Thomas
                    F{var} = [F{var}, l(a)*cross(g{b},g{c}) + l(b)*cross(g{c},g{a}) ...
                        + l(c)*cross(g{a},g{b})];
                    % curls of H(curl) face type 2
                    for i = 0:p-1
                        for j = 0:p-1-i
                            F{var} = [F{var}, cross(grad(vF(j,a,b,c)),grad(uF(i,a,b)))];
                        end
                    end
                    % curls of H(curl) face type 3
                    N0 = g{a}*l(b) - l(a)*g{b};
                    for j = 0:p-1
                        F{var} = [F{var}, rotN(N0*vF(j,a,b,c))];
                    end
                end
            end

            % cell functions (same for both reference tetrahedra)
            uC = @(i) subs(leg.legIntS(i+2),[x,t],[l(1)-l(2),l(1)+l(2)]);
            vC = @(j) l(3)*subs(leg.legS(j+1),[x,t],[2*l(3)-(1-l(4)),1-l(4)]);
            wC = @(k) l(4)*subs(leg.leg(k+1),x,2*l(4)-1);
            N0 = g{1}*l(2) - l(1)*g{2};
            I = sym(zeros(3,0));
            for i = 0:p-2
                for j = 0:p-2-i
                    for k = 0:p-2-i-j
                        I = [I, cross(grad(vC(j)),grad(uC(i)*wC(k))), ...
                            cross(grad(wC(k)),grad(uC(i)*vC(j)))]; %#ok<AGROW>
                    end
                end
            end
            for j = 0:p-2
                for k = 0:p-2-j
                    I = [I, rotN(N0*vC(j)*wC(k))]; %#ok<AGROW>
                end
            end
            for i = 0:p-2
                for j = 0:p-2-i
                    for k = 0:p-2-i-j
                        I = [I, wC(k)*cross(grad(uC(i)),grad(vC(j)))]; %#ok<AGROW>
                    end
                end
            end
            for j = 0:p-2
                for k = 0:p-2-j
                    I = [I, l(2)*wC(k)*cross(g{1},grad(vC(j)))]; %#ok<AGROW>
                end
            end
            for j = 0:p-2
                I = [I, l(2)*vC(j)*cross(g{1},g{4})]; %#ok<AGROW>
            end

            N = cell(1,2);
            divN = cell(1,2);
            for var = 1:2
                B = simplify([F{var}, I]);
                D = sym(zeros(1,size(B,2)));
                for i = 1:size(B,2)
                    D(i) = divergence(B(:,i),dr);
                end
                N{var} = matlabFunction(B,'vars',dr);
                divN{var} = matlabFunction(D,'vars',dr);
            end

            obj.faceDOFs = size(F{1},2)/4;
            obj.interiorDOFs = size(I,2);
            obj.DOFsPerElement = size(F{1},2) + size(I,2);
        end

        function P = evalRef(~,fun,lq,refTet)
            % evaluate a basis cell array at the barycentric point lq for
            % both reference elements and select per element
            c = num2cell(lq);
            P = cat(3,fun{1}(c{:}),fun{2}(c{:}));
            P = P(:,:,refTet);
        end

        function [refTet,dofs,detD,Dk,X] = elementData(obj,phys,pIdx,l)
            g = phys.geometry;
            refTet = g.refTet(pIdx);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = g.detD(:,:,pIdx);
            Dk = g.Dk(:,:,pIdx);
            Nl = obj.dim+1;
            elco = reshape(g.co(:,:,g.el(pIdx,1:Nl)'),[],Nl,numel(pIdx));
            X = cell(1,size(l,2));
            for q = 1:size(l,2)
                % local vertex 1 has barycentric 1-sum(l)
                X{q} = pagemtimes(elco,[1-sum(l(:,q));l(:,q)]);
            end
        end

        function A = toSparse(~,phys,dofs,A)
            Ns = size(A,1);
            I = repmat(dofs,1,Ns)';
            J = repelem(dofs',Ns,1);
            A = sparse(I(:),J(:),A(:),phys.DOFs.Nd,phys.DOFs.Nd);
        end
    end
end

function B = materialTimes(mat,X,q,A)
    % mat * A for scalar/matrix materials or function handles of the
    % coordinates
    if isa(mat,'function_handle')
        mat = double(mat(X{q}));
    end
    if isscalar(mat) || (size(mat,1) == 1 && size(mat,2) == 1)
        B = mat.*A;
    else
        B = pagemtimes(mat,A);
    end
end
