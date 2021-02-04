classdef NedelecIngelReduced < ofem_v2.elements.Finite_Elements & handle
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
		function obj = NedelecIngelReduced(dim, deg)
			%NEDELEC Construct an instance of this class
			%   Detailed explanation goes here
			obj.dim = dim;
			if deg > 3
				warning('ofem_v2:elements:NedelecIngel:Unsupported',...
						'Order only implemented up to 3!\n Setting degree to 3');
				deg = 3;
			end
			obj.degree = deg;
			obj.degreeMass = 3;
			obj.degreeStiff = 3;
			obj.edgeDOFs = 6;
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
					
					E = [];
					F = [];
					I = [];
					
					%% Zeroth order H curl base
					% edge 1 2
					N1 = l1*gradient(l2,dr)-l2*gradient(l1,dr);
					% edge 1 3
					N2 = l1*gradient(l3,dr)-l3*gradient(l1,dr);
					% edge 1 4
					N3 = l1*gradient(l4,dr)-l4*gradient(l1,dr);
					% edge 2 3
					N4 = l2*gradient(l3,dr)-l3*gradient(l2,dr);
					% edge 2 4
					N5 = l2*gradient(l4,dr)-l4*gradient(l2,dr);
					% edge 3 4
					N6 = l3*gradient(l4,dr)-l4*gradient(l3,dr);
					
					E = [N1,N2,N3,N4,N5,N6];
					
					%% First order H curl base also contains grad H1
					if obj.degree > 1
						% face 1 2 3
						N1 = 3*l2*l3*gradient(l1,dr)-gradient(l1*l2*l3,dr);
						N2 = 3*l1*l3*gradient(l2,dr)-gradient(l1*l2*l3,dr);
						% face 1 2 4
						N3 = 3*l2*l4*gradient(l1,dr)-gradient(l1*l2*l4,dr);
						N4 = 3*l1*l4*gradient(l2,dr)-gradient(l1*l2*l4,dr);
						% face 1 3 4
						N5 = 3*l3*l4*gradient(l1,dr)-gradient(l1*l3*l4,dr);
						N6 = 3*l1*l4*gradient(l3,dr)-gradient(l1*l3*l4,dr);
						% face 2 3 4
						N7 = 3*l3*l4*gradient(l2,dr)-gradient(l2*l3*l4,dr);
						N8 = 3*l2*l4*gradient(l3,dr)-gradient(l2*l3*l4,dr);
						F = [N1,N2,N3,N4,N5,N6,N7,N8];
					end
					if obj.degree > 2
						% face 1 2 3
						N1 = 4*l2*l3*(l2-l3)*gradient(l1,dr)-gradient(l1*l2*l3*(l2-l3),dr);
						N2 = 4*l1*l3*(l1-l3)*gradient(l2,dr)-gradient(l1*l2*l3*(l1-l3),dr);
						N3 = 4*l1*l2*(l1-l2)*gradient(l3,dr)-gradient(l1*l2*l3*(l1-l2),dr);
						% face 1 2 4
						N4 = 4*l2*l4*(l2-l4)*gradient(l1,dr)-gradient(l1*l2*l4*(l2-l4),dr);
						N5 = 4*l1*l4*(l1-l4)*gradient(l2,dr)-gradient(l1*l2*l4*(l1-l4),dr);
						N6 = 4*l1*l2*(l1-l2)*gradient(l4,dr)-gradient(l1*l2*l4*(l1-l2),dr);
						%face 1 3 4
						N7 = 4*l3*l4*(l3-l4)*gradient(l1,dr)-gradient(l1*l3*l4*(l3-l4),dr);
						N8 = 4*l1*l4*(l1-l4)*gradient(l3,dr)-gradient(l1*l3*l4*(l1-l4),dr);
						N9 = 4*l1*l3*(l1-l3)*gradient(l4,dr)-gradient(l1*l3*l4*(l1-l3),dr);
						%face 2 3 4
						N10 = 4*l3*l4*(l3-l4)*gradient(l2,dr)-gradient(l2*l3*l4*(l3-l4),dr);
						N11 = 4*l2*l4*(l2-l4)*gradient(l3,dr)-gradient(l2*l3*l4*(l2-l4),dr);
						N12 = 4*l2*l3*(l2-l3)*gradient(l4,dr)-gradient(l2*l3*l4*(l2-l3),dr);
						%interior
						N13 = 4*l2*l3*l4*gradient(l1,dr) - gradient(l1*l2*l3*l4);
						N14 = 4*l1*l3*l4*gradient(l2,dr) - gradient(l1*l2*l3*l4);
						N15 = 4*l1*l2*l4*gradient(l3,dr) - gradient(l1*l2*l3*l4);
						F = [F(:,1:2),N1,N2,N3,F(:,3:4),N4,N5,N6,F(:,5:6),N7,N8,N9,F(:,7:8),N10,N11,N12];
						I = [N13,N14,N15];
					end
					if obj.degree > 3
						warning('ofem_v2:elements:NedelecIngel:Unsupported',...
									'Order only implemented up to 3!\n Setting degree to 3');
						obj.degree = 3;
					end
					N = [E,F,I];
					curlN = N;
					for i = 1:size(curlN,2)
						curlN(:,i) = curl(curlN(:,i),dr);
					end
					
					[~,lc] = obj.quaddata(obj.dim,obj.degreeStiff);
					[~,l] = obj.quaddata(obj.dim,obj.degreeMass);
					N = repmat(N,1,1,size(l,2));
					curlN = repmat(curlN,1,1,size(lc,2));
					for i=1:size(curlN,3)
						curlN(:,:,i) = subs(curlN(:,:,i),dr,lc(:,i)');
					end
					for i=1:size(N,3)
						N(:,:,i) = subs(N(:,:,i),dr,l(:,i)');
					end
					obj.faceDOFs = size(F,2);
					obj.edgeDOFs = size(E,2);
					obj.interiorDOFs = size(I,2);
					N = double(N);
					curlN = double(curlN);
					obj.DOFsPerElement = size(N,2);
					obj.N = N;
					obj.curlN = curlN;					
			end
		end
		
		function [w,l] = quaddata(~,dim,order)
			%QUADDATA returns quadrature points and weights.
			%
			% [w,l]=QUADDATA(dim) returns per row the barycentric coordinates l
			% of the quadrature points and weights w sufficient for the
			% specified element and space dimension given through dim. The
			% quadrature rule is a Gaussian quadrature. The quadarture points l
			% yield barycentric coordinates as they are used by
			% phi and dphi. The weight already contains the jacobian scaling
			% factor for the length/area/volume
			%
			% see also OFEM.NEDELECP.BASIS
			%
			switch dim
				case 2
					w=1/factorial(dim);
					l=repmat(1/(dim+1),1,dim+1);
					coord = [0 1 0; 0 0 1];
					l = coord*l';
				case 3
					switch order
						case 1
							w = 1/6;
							l = [1/4,1/4,1/4]';
						case 2
							w=[1;1;1;1]/24;
							a = 0.1381966011250105151795413165634361;
							b = [a,a,a,1-3*a;a,a,1-3*a,a;a,1-3*a,a,a;1-3*a,a,a,a];
							coord = [0 1 0 0; 0 0 1 0; 0 0 0 1];
							l(:,1) = coord*b(:,1);
							l(:,2) = coord*b(:,2);
							l(:,3) = coord*b(:,3);
							l(:,4) = coord*b(:,4);
						case 3
							w = [-4/5,9/20,9/20,9/20,9/20]/6;
							l = [1/4,1/4,1/4;
								1/6,1/6,1/6;
								1/2,1/6,1/6;
								1/6,1/2,1/6;
								1/6,1/6,1/2]';
					end
			end
			w = w(:);
		end
		
		function S = assembleStiffness(obj,phys,pIdx,A)
			[w,l] = obj.quaddata(obj.dim,obj.degreeStiff);
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
            Ne = size(phys.geometry.el  ,1);
			Nl = size(phys.geometry.el,2);
            
            S=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));

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

            S=S*(1./abs(detD));

            I = repmat(dofs,1,size(S,1))';
            I = I(:);
            J = repelem(dofs',size(S,1),1);
            J = J(:);

            S = sparse(I(:),J(:),S(:),phys.DOFs.Nd,phys.DOFs.Nd);
			figure
            spy(S)
			
		end
		
		function M = assembleMass(obj,phys,pIdx,c)
			[w,l] = obj.quaddata(obj.dim,obj.degreeMass);
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
            Ne = size(phys.geometry.el  ,1);
			Nl = size(phys.geometry.el,2);
            
            M=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));

			if isa(c,'function_handle')
				elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
				for q=1:Nq
					X = elco*(l(q,:)');
					phii = DinvT(:,:,pIdx)*(obj.N(:,:,q).*sign);
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
			figure
            spy(M)
		end
		
		function S = assembleDamping(obj,phys,pIdx,i)
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

