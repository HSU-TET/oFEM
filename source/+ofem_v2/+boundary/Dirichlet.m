classdef Dirichlet < handle
    %   Summary of this class goes here
    %   Detailed explanation goes here
    %% Proteries
    properties
    value; 
    boundary; 
    nodes;
    edges;
    faces;
    interiors;
    index;
    
    u;  
    
    end
    %% Methodes
    methods
        function obj = Dirichlet(value, boundaryName, mesh)
            obj.value = value; 
            
            n = size(mesh.bd, 2);
            for i = 1:n
                if isequal(mesh.bd{1,i},boundaryName)
                    obj.boundary = unique(mesh.bd{2,i});
                    obj.index = i;
                    break;
                end
            end
            obj.nodes = unique(obj.boundary(:));
        end

        function u = loadVector(obj, physicalProblem)
        %dB=dirichlet_loadvector(d,nodes,co)
        %dirichlet returns the Dirichlet-originated part of the load vector.
        
       
        
        %What if there is no M,D, or M_robin?!
            Nco = physicalProblem.geometry.Nco;  % number of nodes
%           dof = physicalProblem.dof;
%           D   = d(co(:,:,nodes));
%           M_robin = d(co(:,:,nodes)); 
      
            physicalProblem.DOFs  = setdiff(physicalProblem.DOFs, obj.boundary);

			obj.u = ones(size(obj.boundary));
			obj.u = sparse(obj.boundary,obj.u,obj.u*obj.value,Nco,1);
                    			
            u = obj.u;
            physicalProblem.DOFs.reduceDOFs(obj);
		             
        end
    end
end
