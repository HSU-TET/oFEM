classdef Dirichlet < handle
    %   Summary of this class goes here
    %   Detailed explanation goes here
    %% Proteries
    properties
    value; 
    boundary; 
    index;
    
	dB; %(S+D+M+M_robin)*boundary  
    
    end
    %% Methodes
    methods
        function obj = Dirichlet(value, boundaryName, mesh)
            obj.value = value; 
            
           n = size(mesh.bd, 2);
            for i = 1:n
                if isequal(mesh.bd{1,i},boundaryName)
                    obj.boundary = unique(mesh.bd{2,i});
                    obj.index = n;
                    break;
                end
            end
        end

        


        
        
        function loadVector(obj, physicalProblem)
        %dB=dirichlet_loadvector(d,nodes,co)
        %dirichlet returns the Dirichlet-originated part of the load vector.
        
       
        
        %What if there is no M,D, or M_robin?!
            Nco = physicalProblem.geometry.Nco;  % number of nodes
%           dof = physicalProblem.dof;
%           D   = d(co(:,:,nodes));
%           M_robin = d(co(:,:,nodes)); 
      
            physicalProblem.DOFs  = setdiff(physicalProblem.DOFs, obj.boundary);

			dB = ones(size(obj.boundary));
			dB = sparse(obj.boundary,dB,dB*obj.value,Nco,1);
                    
                    
            obj.dB= (physicalProblem.S+physicalProblem.D+physicalProblem.M+physicalProblem.M_robin)*dB;
			
			physicalProblem.b = physicalProblem.b - obj.dB;
			physicalProblem.u(obj.boundary) = obj.value;
		 
            
        end
    end
end
