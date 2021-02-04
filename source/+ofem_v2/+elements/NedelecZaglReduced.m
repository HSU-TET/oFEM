classdef NedelecZaglReduced < ofem_v2.elements.Finite_Elements & handle
	%NEDELECZAGL Summary of this class goes here
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
		leg;
    end
	
	methods
		function obj = NedelecZaglReduced(dim, deg)
			%NEDELEC Construct an instance of this class
			%   Detailed explanation goes here
			obj.dim = dim;
			obj.degree = deg;
			obj.degreeMass = 5;
			obj.degreeStiff = 5;
			obj.leg = ofem_v2.tools.LegPoly;
		end
		
		function [N,curlN] = computeBasis(obj)
			switch obj.dim
				case 2
					error('ofem:finiteelement:basis:Unsupported',...
						'2D elements are not implemented yet!');
				case 3
					%% First compute the zeroth order Nedelec elements to construct everything
					syms u v w x t;
					dr = [u,v,w];
					
					%% H1 Basis
					l1 = 1-u-v-w;
					l2 = u;
					l3 = v;
					l4 = w;
					
					l = [l1,l2,l3,l4];
					
					E = [];
					F = [];
					I = [];
					
					%% Zeroth order H curl base
					% edge 1 2
					E{1} = l1*gradient(l2,dr)-l2*gradient(l1,dr);
					% edge 1 3
					E{2} = l1*gradient(l3,dr)-l3*gradient(l1,dr);
					% edge 1 4
					E{3} = l1*gradient(l4,dr)-l4*gradient(l1,dr);
					% edge 2 3
					E{4} = l2*gradient(l3,dr)-l3*gradient(l2,dr);
					% edge 2 4
					E{5} = l2*gradient(l4,dr)-l4*gradient(l2,dr);
					% edge 3 4
					E{6} = l3*gradient(l4,dr)-l4*gradient(l3,dr);
					
					obj.leg.computePolynomials(obj.degree+2);
					F = {};
					F{4} = [];
					I = {};
					
					ke = [2,1;
						 3,1;
						 4,1;
						 3,2;
						 4,2;
						 4,3];
					kf = [3,2,1;
						 4,2,1;
						 4,3,1;
						 4,3,2];
% 					for k = 1:6
% 						for i = 0:obj.degree-1
% 							E{k} = [E{k},gradient(subs(obj.leg.legIntS(i+3),[x,t],[l(ke(k,1))-l(ke(k,2)),l(ke(k,1))+l(ke(k,2))]),dr)];
% 						end
% 					end
					for k = 1:4
						for i = 0:obj.degree-2
							for j = 0:obj.degree-2-i
								m = subs(obj.leg.legIntS(i+3),[x,t],[l(kf(k,1))-l(kf(k,2)),l(kf(k,1))+l(kf(k,2))]);
								n = l(kf(k,3))*subs(obj.leg.legS(j+2),[x,t],[l(kf(k,3))-l(kf(k,1))-l(kf(k,2)),l(kf(k,1))+l(kf(k,2))+l(kf(k,3))]);
% 								F{k} = [F{k},gradient(m*n,dr)];
								F{k} = [F{k},gradient(m,dr)*n-m*gradient(n,dr)];
								if i==0
									F{k} = [F{k},(gradient(l(kf(k,1)),dr)*l(kf(k,2))-l(kf(k,1))*gradient(l(kf(k,2)),dr))*n];
								end
							end
						end
					end
					E = [E{1},E{2},E{3},E{4},E{5},E{6}];
					F = [F{1},F{2},F{3},F{4}];
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
							a = (5-sqrt(5))/20;
							b = (5+3*sqrt(5))/20;
							l = [a,a,a;
								 a,b,a;
								 a,a,b;
								 b,a,a]';
						case 3
							w = [-4/5,9/20,9/20,9/20,9/20]/6;
							l = [1/4,1/4,1/4;
								1/6,1/6,1/6;
								1/2,1/6,1/6;
								1/6,1/2,1/6;
								1/6,1/6,1/2]';
						case 5
							a = (7+sqrt(15))/34;
							b = (7-sqrt(15))/34;
							c = (13+3*sqrt(15))/34;
							d = (13-3*sqrt(15))/34;
							e = (5-sqrt(15))/20;
							f = (5+sqrt(15))/20;
							h = (2665-14*sqrt(15))/226800;
							i = (2665+14*sqrt(15))/226800;
							j = 5/567;
							w = [8/405,h,h,h,h,i,i,i,i,j,j,j,j,j,j];
							l = [1/4,1/4,1/4;
								 a,a,a;
								 a,a,c;
								 a,c,a;
								 c,a,a;
								 b,b,b;
								 b,b,d;
								 b,d,b;
								 d,b,b;
								 e,e,f;
								 e,f,e;
								 f,e,e;
								 e,f,f;
								 f,e,f;
								 f,f,e]';
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

            S=S*ofem_v2.tools.matrixarray((1./abs(detD)));

            I = repmat(dofs,1,size(S,1))';
            I = I(:);
            J = repelem(dofs',size(S,1),1);
            J = J(:);

            S = sparse(I(:),J(:),S(:),phys.DOFs.Nd,phys.DOFs.Nd);
% 			figure
%             spy(S)
			
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
% 			figure
%             spy(M)
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

