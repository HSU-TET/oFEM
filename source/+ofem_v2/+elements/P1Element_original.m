classdef P1Element < ofem_v2.elements.Finite_Elements & handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
  
    
    properties (Access = public)
		shape_function;
        gausianWeight; 
        gausianLength;
        gradShapeFunctions;
        pointwiseProductShapeFunctions;		
    end
    
    properties (Access = public)
        Stiffness_matrix;
        Damping_matrix; 
        Mass_matrix; 
        Load_vector;
    end
    
    methods
        
        function obj = P1Element(mesh)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
			
			obj.calcGausQuad(mesh.dim);     
            obj.calcPointwiseProductShapeFunctions(mesh.dim);
            obj.setShape_function(obj.gausianLength);
            obj.calcGradShapeFunctions(obj.gausianLength);
		end
		
		function assembleRHS(obj, physical_problem, pIdx)        
           
            f = physical_problem.volume_force.value;
            detDLoc= Physical_Problem.geometry.detD(:,:,pIdx);
         
            
            Nl   = size(obj.gausianLength ,2); % number of barycentric coordinates
            Nc   = size(physical_problem.geometry.co ,3);
            Nq   = size(obj.gausianWeight  ,1);
            Ns   = size(obj.shape_function,1);
            Ne   = size(physical_problem.geometry.el(dIdx,:) ,1);

            % elco*l gives global quadrature points => eval there
            elLoc = physical_problem.geometry.el(dIdx,:);
            elco = reshape(physical_problem.geometry.co(:,:,elLoc(:,1:Nl)'),[],Nl,Ne);

            F    = ofem_v2.tools.matrixarray(zeros(1,Ns,Ne));

            for q=1:Nq
                X = elco*(obj.gausianLength(q,:)');
                
                F = F + f(X)*(obj.gausianWeight(q)*obj.shape_function(:,q)'); 
            end

            F  = F*detDLoc;
            elLoc = elLoc';
            obj.Load_vector  = sparse(elLoc(:),1,F(:),Nc,1);
            
        end
        
        function assembleLHS(obj, Physical_Problem, pIdx, pN)
           
           elLoc= Physical_Problem.geometry.el(pIdx,:); 
           co= Physical_Problem.geometry.co;
           detDLoc= Physical_Problem.geometry.detD(:,:,pIdx); % matrixarray(abs(det)) ?
           DinvTLoc= Physical_Problem.geometry.DinvT(:,:,pIdx);
           matpara= 0; 
         
        
           if(Physical_Problem.has_stiffness)
			   if ~isempty(Physical_Problem.geometry.parts{4,pN}.(Physical_Problem.paraS))
                   matpara=Physical_Problem.geometry.parts{4,pN}.(Physical_Problem.paraS);
			   else
				   error('Alles kaputt!!!');
			   end
		   end
            
               
            Ns = size(obj.gradShapeFunctions,1);
            Nq = size(obj.gausianWeight,1);
            Ne = size(elLoc  ,1);
            Nc = size(co  ,3);

           Stiffness_matrix=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
            
            
            if isa(matpara,'function_handle') % function_handle if material is not homogeneous
                Nl   = size(obj.gausianLength,2); % number of barycentric coordinates
                elco = reshape(co(:,:,elLoc(:,1:Nl)'),[],Nl,Ne);

                for q=1:Nq
                    X = elco*(obj.gausianLength(q,:)'); %jeder punkt in globaler sicht
                    globdphi = DinvTLoc*obj.gradShapeFunctions(:,:,q)'; %Achtung abh?ngig von der matrix
                    Stiffness_matrix=Stiffness_matrix+globdphi'*...
                        (obj.gausianWeight(q)* matpara(X))*globdphi;
                end
            else
                for q=1:Nq
                    globdphi = DinvTLoc*obj.gradShapeFunctions(:,:,q)'; %Achtung abh?ngig von der matrix
                    Stiffness_matrix=Stiffness_matrix+globdphi'*...
                        (obj.gausianWeight(q)* matpara)*globdphi;
                end
            end
            
            Stiffness_matrix=Stiffness_matrix*detDLoc;

            J = repmat(1:Ns,Ns,1);
            I = elLoc(:,J')';
            J = elLoc(:,J )';

            obj.Stiffness_matrix{pN} = sparse(I(:),J(:),Stiffness_matrix(:),Nc,Nc); 
            
            if(Physical_Problem.has_damping)
               for i = 1:size(Physical_Problem.geometry.parts,1)
                   para=Physical_Problem.geometry.parts{i,pN};
                   if isequal(para(1,1), {'damping'})
                   matpara= para{1,2}; %TODO
                   end
               end
              
                
            Ns = size(obj.gradShapeFunctions,1);
            Nq = size(obj.gausianWeight,1);
            Ne = size(elLoc  ,1);
            Nc = size(co  ,3);

            obj.Damping_matrix=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
            
            if isa(matpara,'function_handle')
                Nl   = size(obj.gausianLength,2); % number of barycentric coordinates
                elco = reshape(co(:,:,elLoc(:,1:Nl)'),[],Nl,Ne);

                for q=1:Nq
                    X = elco*(obj.gausianLength(q,:)');
                    obj.Damping_matrix=obj.Damping_matrix+obj.shape_functions(:,q)*...
                        ((obj.gausianWeight(q)*matpara(X))'*(DinvTLoc*obj.gradShapeFunctions(:,:,q)'));
                end
            else
                for q=1:Nq
                    obj.Damping_matrix=obj.Damping_matrix+obj.shape_functions(:,q)*...
                        ((obj.gausianWeight(q)*matpara)'*(DinvTLoc*obj.gradShapeFunctions(:,:,q)'));
                end
            end

            obj.Damping_matrix=obj.Damping_matrix*detDLoc;

            J = repmat(1:Ns,Ns,1);
            I = elLoc(:,J')';
            J = elLoc(:,J )';

            obj.Damping_matrix{pN} = sparse(I(:),J(:),obj.Damping_matrix(:),Nc,Nc);
                
            end
            
            if(Physical_Problem.has_mass)
                m = size(Physical_Problem.geometry.parts,1);
               for i = 1:m
                   para=Physical_Problem.geometry.parts{i,pN};
                   if isequal(para(1,1), {'mass'})
                   matpara= para{1,2}; %TODO
                   end
               end
               
            Ns = size(obj.pointwiseProductShapeFunctions,1);
            Nc = size(co  ,3);
      
            Mass_matrix=ofem_v2.tools.matrixarray(zeros(Ns,Ns,Ne));
      
            if isa(matpara,'function_handle') % function_handle if material is not homogeneous
                Nl   = size(obj.gausianLength,2); % number of barycentric coordinates
                elco = reshape(co(:,:,elLoc(:,1:Nl)'),[],Nl,Ne);

                for q=1:Nq
                    X = elco*(obj.gausianLength(q,:)'); %jeder punkt in globaler sicht
                    globdphi = obj.pointwiseProductShapeFunctions(:,:,q)'; 
                    Mass_matrix= Mass_matrix+globdphi'*...
                       (obj.gausianWeight(q)* matpara(X))*globdphi;
                end
            else               
             Mass_matrix= (matpara*obj.pointwiseProductShapeFunctions)*detDLoc;
            end
            J = repmat(1:Ns,Ns,1);
            I = elLoc(:,J')';
            J = elLoc(:,J )';

            obj.Mass_matrix{pN}= sparse(I(:),J(:), Mass_matrix(:),Nc,Nc);
            end
          end      
                
	end	
	methods (Access = protected)
			
		        
        function calcGausQuad(obj,dim)
            obj.gausianWeight =1/factorial(dim);
            obj.gausianLength=repmat(1/(dim+1),1,dim+1);
        end

        function calcPointwiseProductShapeFunctions(obj, dim)
            obj.pointwiseProductShapeFunctions = (ones(dim+1)+eye(dim+1))/factorial(2+dim);
        end
        
        function setShape_function(obj,l)
            obj.shape_function = l';
        end
        
        function calcGradShapeFunctions(obj, l)
            dim = size(l,2)-1;
            obj.gradShapeFunctions = repmat([-1*ones(1,dim);eye(dim)],1,1,size(l,1));
            
        end
      
        function load_vector =volume_force(detD,phi,w,l,f,el,co)
        %volume_force returns the volume force part of the load vector.
        %
        % b=volume_force(detD,phi,w,l,f,el,co) computes the force in
        % terms of a quadrature rule. w and l are the weights and
        % quadrature points of the rule, respectively. phi contains the
        % values of the shape functions evaluated at the quadrature points
        % and f is a functions handle returning the value of the force
        % distribution at arbitrary points.
        %
        % See also ofem.mesh.jacobian_data, ofem.finiteelement.phi,
        % ofem.quassianquadrature.data
        %
            Nl   = size(l  ,2); % number of barycentric coordinates
            Nc   = size(co ,3);
            Nq   = size(w  ,1);
            Ns   = size(phi,1);
            Ne   = size(el ,1);

            % elco*l gives global quadrature points => eval there
            elco = reshape(co(:,:,el(:,1:Nl)'),[],Nl,Ne);

            F    = ofem_v2.tools.matrixarray(zeros(1,Ns,Ne));

            for q=1:Nq
                X = elco*(l(q,:)');
                
                F = F + f(X)*(w(q)*phi(:,q)');
            end

%             F = permute(double(F*detD),[3,2,1]);

            F  = F*detD;
            el = el';
            load_vector  = sparse(el(:),1,F(:),Nc,1);
        end

        function load_vector =pressure(meas,phi,w,l,g,faces,co)
        %pressure returns the pressure-originated part of the load vector.
        %
        % b=pressure(meas,phi,w,l,g,faces,co) computes the pressure
        % originated vector in terms of a quadrature rule. w and l are the
        % weights and quadrature points of the rule, respectively. phi
        % contains the values of the shape functions evaluated at the
        % quadrature points and g is a functions handle returning the
        % pressure at arbitrary points.
        %
        % See also ofem.mesh.jacobian_data, ofem.mesh.neumann,
        % ofem.finiteelement.phi, ofem.quassianquadrature.data
        %
            
            Nl     = size(l  ,2); % number of barycentric coordinates
            Nc     = size(co ,3);
            Nq     = size(w  ,1);
            Ns     = size(phi,1);
            Nf     = size(faces,1);

            % elco*l gives global quadrature points => eval there
            faceco = reshape(co(:,:,faces(:,1:Nl)'),[],Nl,Nf);

            F      = ofem_v2.tools.matrixarray(zeros(1,Ns,Nf));

            for q=1:Nq
                X = faceco*(l(q,:)');
                F = F + g(X)*(w(q)*phi(:,q)');
            end

%             F = permute(double(F*meas),[3,2,1]);

            F     = F*meas;
            faces = faces';
            load_vector     = sparse(faces(:),1,F(:),Nc,1);
        end
    
    end
  
	
end


