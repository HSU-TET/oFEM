classdef (Abstract) Boundary_Condition < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    %% Properties
    
    properties (Abstract)
        value; 
        boundary;
        
    end
    %% Methods
    methods
             
%         function obj = Boundary_Condition(value, boundaryName, mesh)
%             obj.value = value;
%             
%         end
    end
    
    methods (Abstract)
        setBoundaryCondition(obj);
		loadVector(obj,physical);
    end
end

