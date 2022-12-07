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
            obj.shape_function = 'H1';
		end
		
		function [phi,dPhi] = computeBasis(obj)
			switch obj.dim
				case 1
					syms u x t c;
					dr = u;
					
					l1 = 1-u;
					l2 = u;
					
					N = [l1,l2];
					E = [];
					
					leg = ofem_v2.tools.LegPoly;
					leg.computePolynomials(obj.degree+2);
					
					for i = 0:obj.degree-2
						E = [E,subs(leg.legIntS(i+2),[x,t],[N(2)-N(1),N(1)+N(2)])];
					end
					
					phi = [N,E];
					phi = simplify(phi);
					obj.DOFsPerElement = size(phi,2);
					
					dPhi = sym(zeros(2,size(phi,2)));
					for i = 1:size(phi,2)
						dPhi(:,i) = gradient(phi(i),dr);
					end
					
					phiFunc{1} = matlabFunction(phi,'vars',[u]);
					phiFunc{2} = matlabFunction(phi,'vars',[u]);
					dPhiFunc{1} = matlabFunction(dPhi,'vars',[u]);
					dPhiFunc{1} = matlabFunction(dPhi,'vars',[u]);
					obj.nodeDOFs = 1;
					obj.edgeDOFs = obj.degree-1;
					
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
							E = [E,subs(leg.legIntS(i+2),[x,t],[N(ke(k,2))-N(ke(k,1)),N(ke(k,1))+N(ke(k,2))])];
						end
					end
					
					for i = 0:obj.degree-3
						for j = 0:obj.degree-3
							p = subs(leg.legIntS(i+2),[x,t],[N(2)-N(1),N(1)+N(2)]);
							q = N(3)*subs(leg.leg(j+1),2*N(3)-1);
							I = [I,p*q];
						end
					end
					
					phi = [N,E,I];
					phi = simplify(phi);
					obj.DOFsPerElement = size(phi,2);
					
					dPhi = sym(zeros(2,size(phi,2)));
					for i = 1:size(phi,2)
						dPhi(:,i) = gradient(phi(i),dr);
					end
					
					% Quickfix
					phiFunc{1} = matlabFunction(phi,'vars',[u,v]);
					phiFunc{2} = phiFunc{1};
					dPhiFunc{1} = matlabFunction(dPhi,'vars',[u,v]);
					dPhiFunc{2} = dPhiFunc{1};
					obj.nodeDOFs = 1;
					obj.edgeDOFs = obj.degree-1;
					obj.interiorDOFs = 1/2*(obj.degree-2)*(obj.degree-1);
					
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
					
					ke1 = [1,2;1,3;1,4;2,3;2,4;3,4];
					ke2 = [1,2;1,3;1,4;3,2;2,4;3,4];
					
					E1 = [];
					E2 = [];
					F1 = [];
					F2 = [];
					I = [];
					
					for i = 0:obj.degree-2
						for k = 1:6
							E1 = [E1,subs(leg.legIntS(i+2),[x,t],[N(ke1(k,1))-N(ke1(k,2)),N(ke1(k,1))+N(ke1(k,2))])];
							E2 = [E2,subs(leg.legIntS(i+2),[x,t],[N(ke2(k,1))-N(ke2(k,2)),N(ke2(k,1))+N(ke2(k,2))])];
						end
					end

					kf1 = [1,2,3;1,2,4;1,3,4;2,3,4];
					kf2 = [1,3,2;1,2,4;1,3,4;3,2,4];

					for k = 1:4
						for i = 0:obj.degree-3
							for j = 0:obj.degree-i-3
								m = subs(leg.legIntS(i+2),[x,t],[N(kf1(k,1))-N(kf1(k,2)),N(kf1(k,1))+N(kf1(k,2))])*N(kf1(k,3));
								n = subs(leg.legS(j+1),[x,t],[2*N(kf1(k,3))-(N(kf1(k,1))+N(kf1(k,2))+N(kf1(k,3))),N(kf1(k,1))+N(kf1(k,2))+N(kf1(k,3))]);
								F1 = [F1,m*n];
								m = subs(leg.legIntS(i+2),[x,t],[N(kf2(k,1))-N(kf2(k,2)),N(kf2(k,1))+N(kf2(k,2))])*N(kf2(k,3));
								n = subs(leg.legS(j+1),[x,t],[2*N(kf2(k,3))-(N(kf2(k,1))+N(kf2(k,2))+N(kf2(k,3))),N(kf2(k,1))+N(kf2(k,2))+N(kf2(k,3))]);
								F2 = [F2,m*n];
							end
						end
					end

					for i = 0:obj.degree-4
						for j = 0:obj.degree-i-4
							for k = 0:obj.degree-i-j-4
								m = subs(leg.legIntS(i+2),[x,t],[N(1)-N(2),N(1)+N(2)])*N(3);
								n = subs(leg.legS(j+1),[x,t],[2*N(3)-(1-N(4)),1-N(4)])*N(4);
								p = subs(leg.leg(k+1),[x],[2*N(4)-1]);
								I = [I,m*n*p];
							end
						end
					end
					
					N = reshape([N,N],1,size(N,2),2);
                    E = reshape([E1,E2],1,size(E1,2),2);
                    F = reshape([F1,F2],1,size(F1,2),2);
                    I = reshape([I,I],1,size(I,2),2);
					
					phi = [N,E,F,I];
					phi = simplify(phi);
					obj.DOFsPerElement = size(phi,2);
					
					dPhi = sym(zeros(3,size(phi,2)));
					for i = 1:size(phi,2)
						dPhi(:,i,1) = gradient(phi(:,i,1),dr);
						dPhi(:,i,2) = gradient(phi(:,i,2),dr);
					end
					
					phiFunc{1} = matlabFunction(phi(:,:,1),'vars',[u,v,w]);
					phiFunc{2} = matlabFunction(phi(:,:,2),'vars',[u,v,w]);
					dPhiFunc{1} = matlabFunction(dPhi(:,:,1),'vars',[u,v,w]);
					dPhiFunc{2} = matlabFunction(dPhi(:,:,2),'vars',[u,v,w]);
					
					obj.nodeDOFs = 1;
					obj.edgeDOFs = size(E,2)/6;
					obj.faceDOFs = size(F,2)/4;
					obj.interiorDOFs = size(I,2);
					
			end
			
			
			phi = phiFunc;
			dPhi = dPhiFunc;
			obj.phi = phi;
			obj.dPhi = dPhi;
		end
		
		function S = assembleStiffness(obj, phys, pIdx, mat)
			refTet = phys.geometry.refTet(pIdx);
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeStiff);
			dofs = phys.DOFs.el2DOF(pIdx,:);
			detD = phys.geometry.detD(:,:,pIdx);
			DinvT = phys.geometry.DinvT(:,:,pIdx);
			
			Ns = obj.DOFsPerElement;
			Nq = size(w   ,1);
			Ne = length(pIdx);
			Nl = size(phys.geometry.el,2);
			
			S = ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));

			chver = split(version,'.');
			chver = str2double(chver{2});
            if chver >= 9
                S = double(S);
                detD = double(detD);
                DinvT = double(DinvT);
            end
			
			if isa(mat,'function_handle')
				elco = reshape(phys.geometry.co(:,:,phys.geometry.el(pIdx,1:Nl)'),[],Nl,Ne);
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					dphi(:,:,1) = obj.dPhi{1}(lTemp{:});
					dphi(:,:,2) = obj.dPhi{2}(lTemp{:});
					X = elco*[l(:,q);1-sum(l(:,q))];
                    if chver < 9
					    dphi = DinvT*ofem_v2.tools.matrixarray(dphi(:,:,refTet));
					    S = S+w(q)*(dphi'*mat(X)*dphi);
                    else
                        dphi = pagemtimes(DinvT,dphi(:,:,refTet));
					    S = S+w(q)*pagemtimes(dphi,'transpose',pagemtimes(mat(X),dphi),'none');
                    end
				end
			else
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					dphi(:,:,1) = obj.dPhi{1}(lTemp{:});
					dphi(:,:,2) = obj.dPhi{2}(lTemp{:});
					if chver < 9
                        dphi =  DinvT*ofem_v2.tools.matrixarray(dphi(:,:,refTet));
                        S = S+w(q)*(dphi'*mat*dphi);
                    else
                        dphi = pagemtimes(DinvT,dphi(:,:,refTet));
                        S = S+w(q)*pagemtimes(dphi,'transpose',mat*dphi,'none');
					end
				end
			end
			
			if chver < 9
				S=S*ofem_v2.tools.matrixarray(abs(detD));
			else
				S = pagemtimes(S,abs(detD));
			end
			
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
			
			Ns = obj.DOFsPerElement;
			Nq = size(w   ,1);
			Ne = length(pIdx);
			Nl = size(phys.geometry.el,2);
			
			M=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
			
			chver = split(version,'.');
			chver = str2double(chver{2});
            if chver >= 9
                M = double(M);
                detD = double(detD);
            end

			if isa(mat,'function_handle')
				elco = reshape(phys.geometry.co(:,:,phys.geometry.el(pIdx,1:Nl)'),[],Nl,Ne);
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi(:,:,1) = obj.phi{1}(lTemp{:});
					phi(:,:,2) = obj.phi{2}(lTemp{:});
					X = elco*[l(:,q);1-sum(l(:,q))];
                    if chver < 9
					    phi = ofem_v2.tools.matrixarray(phi(:,:,refTet));
					    M = M+w(q)*(phi'*mat(X)*phi);
                    else
                        phi = phi(:,:,refTet);
					    M = M+w(q)*pagemtimes(phi,'transpose',pagemtimes(mat(X),phi),'none');
                    end
				end
			else
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi(:,:,1) = obj.phi{1}(lTemp{:});
					phi(:,:,2) = obj.phi{2}(lTemp{:});
					if chver < 9
						phi = ofem_v2.tools.matrixarray(phi);
                        M = M+w(q)*(phi(:,:,refTet)'*mat*phi(:,:,refTet));
                    else
                        M = M+w(q)*pagemtimes(phi(:,:,refTet),'transpose',mat*phi(:,:,refTet),'none');
					end
				end
			end
			
			if chver < 9
				M = M*ofem_v2.tools.matrixarray(abs(detD));
			else
				M = pagemtimes(M,abs(detD));
			end
			
			I = repmat(dofs,1,size(M,1))';
			%I = I(:);
			J = repelem(dofs',size(M,1),1);
			%J = J(:);
			
			M = sparse(I(:),J(:),M(:),phys.DOFs.Nd,phys.DOFs.Nd);
		end
		
		function D = assembleDamping(obj,phys,pIdx,mat)
		end
		
		function b = volumeForce(obj,phys,value,pIdx)
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim,obj.degreeMass);
			dofs = phys.DOFs.el2DOF(pIdx,:);
			detD = phys.geometry.detD(:,:,pIdx);
			
			Ns = obj.DOFsPerElement;
			Nq = size(w   ,1);
			Ne = length(pIdx);
			Nl = size(phys.geometry.el,2);
			
			F = ofem_v2.tools.matrixarray(zeros(Ns,1,Ne));
			
			
			if isa(value,'function_handle')
				loco = reshape(phys.geometry.co(:,:,phys.geometry.el(pIdx,:)'),[],Nl,Ne);
				for q=1:Nq
					X = loco*([l(:,q);1-sum(l(:,q))]);
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi = obj.phi{1}(lTemp{:});
					F = F+w(q)*(phi'*value(X));
				end
			elseif isa(value,'ofem_v2.tools.matrixarray')
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi = obj.phi{1}(lTemp{:});
					F = F+w(q)*(phi'*value);
				end
			else
				for q=1:Nq
					cnt = ones(size(l(:,q),1),1);
					lTemp = mat2cell(l(:,q),cnt);
					phi = obj.phi{1}(lTemp{:});
					F = F+w(q)*(phi'*value);
				end
			end
			
			F = F*abs(detD);
			
			b = sparse(dofs',1,F(:),phys.DOFs.Nd,1);
		end
	end
end






























