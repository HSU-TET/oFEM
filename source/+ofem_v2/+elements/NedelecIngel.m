classdef NedelecIngel < ofem_v2.elements.Finite_Elements & handle
	%NEDELEC Summary of this class goes here
	%   Detailed explanation goes here
	
    properties (Access = public)
		shape_function;
        gausianWeight; 
        gausianLength;
        gradShapeFunctions;
        pointwiseProductShapeFunctions;		
    end
    
    properties (Access = public)
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
    end
	
	methods
		function obj = NedelecIngel(dim, deg)
			%NEDELEC Construct an instance of this class
			%   Detailed explanation goes here
			obj.dim = dim;
			if deg > 3
				warning('ofem_v2:elements:NedelecIngel:Unsupported',...
						'Order only implemented up to 3!\n Setting degree to 3');
				deg = 3;
			end
			obj.degree = deg;
            deg = 2*(deg+2);
            %if deg == 0; deg = 1; end
			obj.degreeMass = deg;
			obj.degreeStiff = deg;
			obj.edgeDOFs = 6*deg;
			obj.faceDOFs = 4*deg*(deg-1);
			obj.interiorDOFs = deg*(deg-1)*(deg-2)/2;
			obj.DOFsPerElement = obj.edgeDOFs+obj.faceDOFs+obj.interiorDOFs;
		end
		
		function [N,curlN] = computeBasis(obj)
			switch obj.dim
				case 2
					error('ofem:finiteelement:basis:Unsupported',...
						'2D elements are not implemented yet!');
				case 3
					%% First compute the zeroth order Nedelec elements to construct everything
					syms u v w ;
					dr = [u,v,w];
					
					%% H1 Basis
					l1 = 1-u-v-w;
					l2 = u;
					l3 = v;
					l4 = w;
                    
                    l(1) = l1;
                    l(2) = l2;
                    l(3) = l3;
                    l(4) = l4;
					
					E = [];
					F = [];
					I = [];
					
					%% Zeroth order H curl base
					% edge 1 2
                    ke = [1,2;2,3;3,1;2,4;4,3;1,4];
					E{1} = l(ke(1,1))*gradient(l(ke(1,2)),dr)-l(ke(1,2))*gradient(l(ke(1,1)),dr);
					% edge 1 3
					E{2} = l(ke(2,1))*gradient(l(ke(2,2)),dr)-l(ke(2,2))*gradient(l(ke(2,1)),dr);
					% edge 1 4
					E{3} = l(ke(3,1))*gradient(l(ke(3,2)),dr)-l(ke(3,2))*gradient(l(ke(3,1)),dr);
					% edge 2 3
					E{4} = l(ke(4,1))*gradient(l(ke(4,2)),dr)-l(ke(4,2))*gradient(l(ke(4,1)),dr);
					% edge 2 4
					E{5} = l(ke(5,1))*gradient(l(ke(5,2)),dr)-l(ke(5,2))*gradient(l(ke(5,1)),dr);
					% edge 3 4
					E{6} = l(ke(6,1))*gradient(l(ke(6,2)),dr)-l(ke(6,2))*gradient(l(ke(6,1)),dr);
					                    
                    %ke = [1,2;2,3;3,1;2,4;4,3;1,4];
                    
                    for k = 1:6
                        for i = 1:obj.degree
                            E{k} = [E{k},gradient(l(ke(k,1))*l(ke(k,2)),dr)];
                            if obj.degree == 2
                                E{k} = [E{k},gradient(l(ke(k,1))*l(ke(k,2))*(l(ke(k,1))-l(ke(k,2))),dr)];
                            end
                            if obj.degree == 3
                                E{k} = [E{k},gradient(l(ke(k,1))*l(ke(k,2))*(l(ke(k,1))^2-3*l(ke(k,1))*l(ke(k,2))+l(ke(k,2))^2),dr)];
                            end
                        end
                    end
                    
                    F = {};
                    F{4} = [];
                    
                    kf = [1,2,3;1,3,4;1,4,2;2,4,3];
                    
                    for k = 1:4
                        for i = 1:obj.degree
                            F{k} = [F{k},3*l(kf(k,2))*l(kf(k,3))*gradient(l(kf(k,1)),dr)-gradient(l(kf(k,1))*l(kf(k,2))*l(kf(k,3)),dr)];
                            F{k} = [F{k},3*l(kf(k,3))*l(kf(k,1))*gradient(l(kf(k,2)),dr)-gradient(l(kf(k,1))*l(kf(k,2))*l(kf(k,3)),dr)];
                        end
                    end
                    
                    E = [E{1},E{2},E{3},E{4},E{5},E{6}];
                    F = [F{1},F{2},F{3},F{4}];
                    obj.edgeDOFs = size(E,2);
                    obj.faceDOFs = size(F,2);
                    obj.interiorDOFs = size(I,2);
					N = [E,F,I];
					curlN = N;
					for i = 1:size(curlN,2)
						curlN(:,i) = curl(curlN(:,i),dr);
					end
					
					[~,lc] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeStiff);
					[~,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
					N = repmat(N,1,1,size(l,2));
					curlN = repmat(curlN,1,1,size(lc,2));
					for i=1:size(curlN,3)
						curlN(:,:,i) = subs(curlN(:,:,i),dr,lc(:,i)');
					end
					for i=1:size(N,3)
						N(:,:,i) = subs(N(:,:,i),dr,l(:,i)');
					end
					N = double(N);
					curlN = double(curlN);
					obj.DOFsPerElement = size(N,2);
					obj.N = N;
					obj.curlN = curlN;
			end
        end
		
		function S = assembleStiffness(obj,phys,pIdx,A)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeStiff);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            edsign = phys.geometry.el2edsign(pIdx,:);
            fasign = phys.geometry.el2fasign(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            Dk = phys.geometry.Dk(:,:,pIdx);
            
            sign = [repelem(edsign,3,phys.element.edgeDOFs/6,1),...
                repelem(fasign,3,phys.element.faceDOFs/4,1),...
                ones(size(edsign,1)*3,phys.element.interiorDOFs)];
            sign = ofem_v2.tools.matrixarray(sign);
            sign = reshape(sign',phys.element.DOFsPerElement,3,[])';
            
            Ns = size(obj.curlN,2);
            Nq = size(w   ,1);
            Ne = length(pIdx);
            Nl = size(phys.geometry.el,2);
            
            S = ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
            
            if isa(A,'function_handle')
                elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
                for q=1:Nq
                    X = elco*(l(q,:)');
                    dphii = Dk(:,:,pIdx)*(obj.curlN(:,:,q).*sign);
                    S = S+w(q)*(dphii'*A(X)*dphii);
                end
            else
                for q=1:Nq
                    dphii = Dk*(obj.curlN(:,:,q).*sign);
                    S = S+w(q)*(dphii'*A*dphii);
                end
            end
            
            S=S*ofem_v2.tools.matrixarray((1./abs(detD)));
            
            I = repmat(dofs,1,size(S,1))';
            I = I(:);
            J = repelem(dofs',size(S,1),1);
            J = J(:);
            
            S = sparse(I(:),J(:),S(:),phys.DOFs.Nd,phys.DOFs.Nd);
			
		end
		
		function M = assembleMass(obj,phys,pIdx,c)
            [w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
            dofs = phys.DOFs.el2DOF(pIdx,:);
            edsign = phys.geometry.el2edsign(pIdx,:);
            fasign = phys.geometry.el2fasign(pIdx,:);
            detD = phys.geometry.detD(:,:,pIdx);
            DinvT = phys.geometry.DinvT(:,:,pIdx);
            
            sign = [repelem(edsign,3,phys.element.edgeDOFs/6,1),...
                repelem(fasign,3,phys.element.faceDOFs/4,1),...
                ones(size(edsign,1)*3,phys.element.interiorDOFs)];
            sign = ofem_v2.tools.matrixarray(sign);
            sign = reshape(sign',phys.element.DOFsPerElement,3,[])';
            
            Ns = size(obj.curlN,2);
            Nq = size(w   ,1);
            Ne = length(pIdx);
            Nl = size(phys.geometry.el,2);
            
            M=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
            
            if isa(c,'function_handle')
                elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
                for q=1:Nq
                    X = elco*(l(q,:)');
                    phii = DinvT(:,:,pIdx)'*(obj.N(:,:,q).*sign);
                    M = M+w(q)*(phii'*c(X)*phii);
                end
            else
                for q=1:Nq
                    phii = DinvT*(obj.N(:,:,q).*sign);
                    M = M+w(q)*(phii'*c*phii);
                end
            end
            
            M=M*abs(detD);
            
            I = repmat(dofs,1,size(M,1))';
            I = I(:);
            J = repelem(dofs',size(M,1),1);
            J = J(:);
            
            M = sparse(I(:),J(:),M(:),phys.DOFs.Nd,phys.DOFs.Nd);
            % 			figure
            %             spy(M)
		end
		
		function S = assembleDamping(obj,phys,pIdx,i)
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
            
            if obj.degree >= 1
                sten1 = sparse(diag([0,0]));
                %sten2 = sparse(diag([1,0]));
                
                %eye1 = speye(obj.degree);
                %eye2 = speye(sum(1:obj.degree-2));
                
                %stencil = blkdiag(kron(eye1,sten1),kron(eye2,sten2));
                
                stencil = sten1;
                
                Y_F = kron(speye(Nfa),stencil);
            else
                Y_F = [];
            end
            
            if obj.degree >= 2
                sten1 = sparse(diag([1,0,0,0]));
                sten2 = sparse(diag([1,0,0]));
                
                eye1 = speye(obj.degree-2);
                eye2 = speye(sum(1:obj.degree-3));
                
                stencil = blkdiag(kron(eye1,sten1),kron(eye2,sten2));
                
                Y_I = kron(speye(Nint),stencil);
            else
                Y_I = [];
            end
            
            Y_g = blkdiag(Y_N0,Y_E,Y_F,Y_I);
            S_g = Y_g'*phys.M*Y_g;
        end
		
	end
	
	methods (Access = protected)
		function calcGradShapeFunctions(shape_function)
			test=1;
		end
		function calcPointwiseProductShapeFunctions(shape_function)
			test = 1;
		end
	end
end

