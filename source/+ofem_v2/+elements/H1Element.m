classdef H1Element < ofem_v2.elements.Finite_Elements & handle
	%UNTITLED5 Summary of this class goes here
	%   Detailed explanation goes here
	
	
	properties (Access = public)
		shape_function;
		phi;
		dPhi;
		dim;
		degree;
		Stiffness_matrix;
		Damping_matrix;
		Mass_matrix;
		Load_vector;
		nodeDOFs = 0;
		edgeDOFs = 0;
		faceDOFs = 0;
		interiorDOFs = 0;
		DOFsPerElement = 0;
		degreeMass;
		degreeStiff;
	end
	
	methods
		
		function obj = H1Element(dim,deg)
			%UNTITLED5 Construct an instance of this class
			%   Detailed explanation goes here
			
			obj.dim = dim;
			obj.degree = deg;
			obj.degreeMass = deg;
			obj.degreeStiff = deg;
		end
		
		function [phi,dPhi] = computeBasis(obj)
			switch obj.dim
				case 2
					syms u v x t c;
					dr = [u,v];
					
					l1 = 1-u-v;
					l2 = u;
					l3 = v;
					
					N = [l1,l2,l3];
					E = [];
					I = [];
					
					ke = [1,2;1,3;2,3];
					
					leg = ofem_v2.tools.LegPoly;
					leg.computePolynomials(obj.degree+2);
					
					for i = 0:obj.degree-2
						for k = 1:3
							E = [E,subs(leg.legIntS(i+2),[x,t],[N(ke(k,1))-N(ke(k,2)),N(ke(k,1))+N(ke(k,2))])];
						end
					end
					
					phi = [N,E];
					phi = simplify(phi);
					obj.DOFsPerElement = size(phi,2);
					
					dPhi = sym(zeros(2,size(phi,2)));
					for i = 1:size(phi,2)
						dPhi(:,i) = gradient(phi(i),dr);
					end
					
					phiFunc = matlabFunction(phi,'vars',[u,v]);
					dPhiFunc = matlabFunction(dPhi,'vars',[u,v]);
					obj.nodeDOFs = 1;
					obj.edgeDOFs = obj.dim-1;
					obj.interiorDOFs = 1/2*(obj.dim-2)*(obj.dim-1);
					
				case 3
					syms u v w x t c;
					dr = [u,v,w];
					
					l1 = 1-u-v-w;
					l2 = u;
					l3 = v;
					l4 = w;
					
					N = [l1,l2,l3,l4];
					
					leg = ofem_v2.tools.LegPoly;
					leg.computePolynomials(obj.degree+2);
					
					ke = [1,2;1,3;1,4;2,3;2,4;3,4];
					
					E = [];
					F = [];
					I = [];
					
					for i = 0:obj.degree-1
						for k = 1:6
							E = [E,subs(leg.legIntS(i+2),[x,t],[N(ke(k,1))-N(ke(k,2)),N(ke(k,1))+N(ke(k,2))])];
						end
					end
					
					phi = [N,E];
					phi = simplify(phi);
					obj.DOFsPerElement = size(phi,2);
					
					dPhi = sym(zeros(3,size(phi,2)));
					for i = 1:size(phi,2)
						dPhi(:,i) = gradient(phi(i),dr);
					end
					
					phiFunc = matlabFunction(phi,'vars',[u,v,w]);
					dPhiFunc = matlabFunction(dPhi,'vars',[u,v,w]);
					
					obj.nodeDOFs = 1;
					obj.edgeDOFs = obj.dim-1;
					obj.faceDOFs = 1/2*(obj.dim-2)*(obj.dim-1);
					obj.interiorDOFs = 1/6*(obj.dim-3)*(obj.dim-2)*(obj.dim-1);
					
			end

			
			phi = phiFunc;
			dPhi = dPhiFunc;
			obj.phi = phi;
			obj.dPhi = dPhi;
		end
		
		function S = assembleStiffness(obj, phys, pIdx, mat)
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeStiff);
			dofs = phys.DOFs.el2DOF(pIdx,:);
			detD = phys.geometry.detD(:,:,pIdx);
			DinvT = phys.geometry.DinvT(:,:,pIdx);
			
			Ns = obj.DOFsPerElement;
			Nq = size(w   ,1);
			Ne = length(pIdx);
			Nl = size(phys.geometry.el,2);
			
			S = ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
			
			if isa(mat,'function_handle')
				elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
				for q=1:Nq
					X = elco*(l(q,:)');
					dphi = DinvT(:,:,pIdx)*obj.dPhi(l(1,q),l(2,q),l(3,q));
					S = S+w(q)*(dphi'*mat(X)*dphi);
				end
			else
				for q=1:Nq
					%dphi = DinvT*obj.dPhi(l(1,q),l(2,q),l(3,q));
					lTemp = mat2cell(l(:,q),[1,1]);
					dphi = DinvT*obj.dPhi(lTemp{:});
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
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
			dofs = phys.DOFs.el2DOF(pIdx,:);
			detD = phys.geometry.detD(:,:,pIdx);
			
			Ns = obj.DOFsPerElement;
			Nq = size(w   ,1);
			Ne = length(pIdx);
			Nl = size(phys.geometry.el,2);
			
			M=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
			
			if isa(mat,'function_handle')
				elco = reshape(phys.geometry.co(:,:,el(pIdx,1:Nl)'),[],Nl,Ne);
				for q=1:Nq
					X = elco*(l(q,:)');
					phi = obj.phi{1}(l(1,q),l(2,q),l(3,q));
					M = M+w(q)*(phi'*mat(X)*phi);
				end
			else
				for q=1:Nq
					phi = obj.phi{1}(l(1,q),l(2,q),l(3,q));
					M = M+w(q)*(phi'*mat*phi);
				end
			end
			
			M=M*abs(detD);
			
			I = repmat(dofs,1,size(M,1))';
			%I = I(:);
			J = repelem(dofs',size(M,1),1);
			%J = J(:);
			
			M = sparse(I(:),J(:),M(:),phys.DOFs.Nd,phys.DOFs.Nd);
		end
		
		function D = assembleDamping(obj,phys,pIdx,mat)
		end
		
	end
end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
