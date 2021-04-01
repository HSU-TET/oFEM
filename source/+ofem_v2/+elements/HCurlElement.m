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
        edgeDOFs;
        faceDOFs;
        interiorDOFs;
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
            obj.degreeMass = 2*(deg);
            obj.degreeStiff = 2*(deg);
        end
        
        function [N,curlN] = computeBasis(obj)
            switch obj.dim
                case 2
                    error('ofem:finiteelement:basis:Unsupported',...
                        '2D elements are not implemented yet!');
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
                    
                    E = {};
                    
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
            refTet = phys.geometry.refTet(pIdx);
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeStiff);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            Dk = phys.geometry.Dk(:,:,pIdx);
            
            Ns = obj.DOFsPerElement;
            Nq = size(w   ,1);
            Ne = length(pIdx);
            Nl = size(phys.geometry.el,2);
            
            S = ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
            
            if isa(mat,'function_handle')
                elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
                for q=1:Nq
                    X = elco*(l(q,:)');
                    dphii = Dk(:,:,pIdx)*(obj.curlN(:,:,q).*sign);
                    S = S+w(q)*(dphii'*mat(X)*dphii);
                end
            else
                for q=1:Nq
                    dphi(:,:,1) = obj.curlN{1}(l(1,q),l(2,q),l(3,q));
                    dphi(:,:,2) = obj.curlN{2}(l(1,q),l(2,q),l(3,q));
                    dphi =  Dk*ofem_v2.tools.matrixarray(dphi(:,:,refTet));
%                     dphii = (Dk*obj.curlN(:,:,q));
                    S = S+w(q)*(dphi'*mat*dphi);
                end
            end
            
            S=S*ofem_v2.tools.matrixarray(1./abs(detD));
            
            I = repmat(dofs,1,size(S,1))';
            I = I(:);
            J = repelem(dofs',size(S,1),1);
            J = J(:);
            
            S = sparse(I(:),J(:),S(:),phys.DOFs.Nd,phys.DOFs.Nd);
            
        end
        
        function M = assembleMass(obj,phys,pIdx,mat)
            refTet = phys.geometry.refTet(pIdx);
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            DinvT = phys.geometry.DinvT(:,:,pIdx);
                        
            Ns = obj.DOFsPerElement;
            Nq = size(w   ,1);
            Ne = length(pIdx);
            Nl = size(phys.geometry.el,2);
            
            M=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
            
            if isa(mat,'function_handle')
                elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
                for q=1:Nq
                    X = elco*(l(q,:)');
                    phii = DinvT(:,:,pIdx)'*(obj.N(:,:,q).*sign);
                    M = M+w(q)*(phii'*mat(X)*phii);
                end
            else
                for q=1:Nq
                    phi(:,:,1) = obj.N{1}(l(1,q),l(2,q),l(3,q));
                    phi(:,:,2) = obj.N{2}(l(1,q),l(2,q),l(3,q));
                    phi =  DinvT*ofem_v2.tools.matrixarray(phi(:,:,refTet));
                    M = M+w(q)*(phi'*mat*phi);
                end
            end
            
            M=M*abs(detD);
            
            I = repmat(dofs,1,size(M,1))';
            %I = I(:);
            J = repelem(dofs',size(M,1),1);
            %J = J(:);
            
            M = sparse(I(:),J(:),M(:),phys.DOFs.Nd,phys.DOFs.Nd);
            
            % 			figure
            %             spy(M)
        end
        
        function D = assembleDamping(obj,phys,pIdx,mat)
        end
        
        function D = assembleVelocityTerm(obj,phys,pIdx,mat,v)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            DinvT = phys.geometry.DinvT(:,:,pIdx);
            Dk = phys.geometry.Dk(:,:,pIdx);
            
            Ns = obj.DOFsPerElement;
            Nq = size(w   ,1);
            Ne = length(pIdx);
            Nl = size(phys.geometry.el,2);
            
            D = ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
            
            refTet = phys.geometry.refTet(pIdx);
            
            if isa(mat,'function_handle')
                elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
                for q=1:Nq
                    X = elco*(l(q,:)');
                    phii = DinvT(:,:,pIdx)'*(obj.N(:,:,q).*sign);
                    M = M+w(q)*(phii'*mat(X)*phii);
                end
            else
                for q=1:Nq
                    phi(:,:,1) = obj.N{1}(l(1,q),l(2,q),l(3,q));
                    phi(:,:,2) = obj.N{2}(l(1,q),l(2,q),l(3,q));
                    dphi(:,:,1) = obj.curlN{1}(l(1,q),l(2,q),l(3,q));
                    dphi(:,:,2) = obj.curlN{2}(l(1,q),l(2,q),l(3,q));
                    phi =  DinvT*ofem_v2.tools.matrixarray(phi(:,:,refTet));
                    dphi =  cross(repmat(v,1,Ns,length(pIdx)),Dk*ofem_v2.tools.matrixarray(dphi(:,:,refTet)));
                    D = D+w(q)*(dphi'*mat*phi);
                end
            end
            
            I = repmat(dofs,1,size(D,1))';
            I = I(:);
            J = repelem(dofs',size(D,1),1);
            J = J(:);
            
            D = sparse(I(:),J(:),D(:),phys.DOFs.Nd,phys.DOFs.Nd);
        end
                
        function F = volumeForce(obj,phys,f,pIdx)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            DinvT = phys.geometry.DinvT(:,:,pIdx);
            
            Ns = obj.DOFsPerElement;
            Nq = size(w   ,1);
            Ne = length(pIdx);
            Nl = size(phys.geometry.el,2);
            
            F = ofem_v2.tools.matrixarray(zeros(1,Ns,Ne));
            
            refTet = phys.geometry.refTet(pIdx);
            
            if isa(f,'function_handle')
                elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
                for q=1:Nq
                    X = elco*(l(q,:)');
                    A = f(X);
                    phii = DinvT(:,:,pIdx)*(obj.N(:,:,q).*sign);
                    F = F+w(q)*(A'*phii);
                end
            elseif isa(f,'ofem_v2.tools.matrixarray')
                for q=1:Nq
                    phi(:,:,1) = obj.N{1}(l(1,q),l(2,q),l(3,q));
                    phi(:,:,2) = obj.N{2}(l(1,q),l(2,q),l(3,q));
                    phi =  DinvT*ofem_v2.tools.matrixarray(phi(:,:,refTet));
                    F = F+w(q)*(f'*phi);
%                 % accepts rhs data converted to integration points of the
%                 % shape dim(f) dim(q) dim(el)
%                     phi(:,:,1) = obj.curlN{1}(l(1,q),l(2,q),l(3,q));
%                     phi(:,:,2) = obj.curlN{2}(l(1,q),l(2,q),l(3,q));
%                     phi =  DinvT*ofem_v2.tools.matrixarray(phi(:,:,refTet));
%                     F = F+w(q)*(f(:,q,:)'*phi);
                end
            else
                f = ofem_v2.tools.matrixarray(repmat(f,1,1,Ne));
                for q=1:Nq
                    phi(:,:,1) = obj.N{1}(l(1,q),l(2,q),l(3,q));
                    phi(:,:,2) = obj.N{2}(l(1,q),l(2,q),l(3,q));
                    phi =  DinvT*ofem_v2.tools.matrixarray(phi(:,:,refTet));
                    F = F+w(q)*(f'*phi);
                end
            end
            
            F = F'*abs(detD);
            
            I = dofs';
            I = I(:);
            
            F = sparse(I(:),1,F(:),phys.DOFs.Nd,1);
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
            l(4,:) = 1-sum(l,1);
            Nl = size(l,1);
            Np = size(l,2);
            elco = reshape(phys.geometry.co(:,:,phys.geometry.el(:,:)'),[],Nl,Ne);
            X = zeros(3,Ne,Np);
            for i = 1:Np
                X(:,:,i) = reshape(elco*l(:,i),3,[]);
            end
            elu = ofem_v2.tools.matrixarray(reshape(u(phys.DOFs.el2DOF(:,:)'),[],1,Ne));
            for q = 1:Np
                dphi(:,:,1) = phys.element.curlN{1}(l(1,q),l(2,q),l(3,q));
                dphi(:,:,2) = phys.element.curlN{2}(l(1,q),l(2,q),l(3,q));
                dphi = ofem_v2.tools.matrixarray(dphi(:,:,phys.geometry.refTet));
                U(:,:,q) = reshape((Dk*dphi)*elu*ofem_v2.tools.matrixarray(1/(detD)),3,[]);
            end
        end

        function [S_g,Y_g] = AMS(obj,phys)
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
            
            if obj.degree > 1 && obj.faceDOFs>0
                sten1 = sparse(diag([1,0,0]));
                sten2 = sparse(diag([1,0]));
                
                eye1 = speye(obj.degree-1);
                eye2 = speye(sum(1:obj.degree-2));
                
                stencil = blkdiag(kron(eye1,sten1),kron(eye2,sten2));
                                
                Y_F = kron(speye(Nfa),stencil);
            else
                Y_F = [];
            end
            
            if obj.degree > 2 && obj.interiorDOFs > 0
                sten1 = sparse(diag([1,0,0,0]));
                sten2 = sparse(diag([1,0,0]));
                
                eye1 = speye(sum(1:obj.degree-2));
                eye2 = speye(sum(sum(triu(repelem(1:obj.degree-3,obj.degree-3,1)'))));
                
                stencil = blkdiag(kron(eye1,sten1),kron(eye2,sten2));
                
                Y_I = kron(speye(Nint),stencil);
            else
                Y_I = [];
            end
            
            Y_g = blkdiag(Y_N0,Y_E,Y_F,Y_I);
            S_g = Y_g'*phys.M*Y_g;
		end
        
    end
    
end

