classdef Initial_Condition < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    % Wenn der solver f�r parabolic anfangen will, setzt er die �brigen
    % dofs erstmal auf initial_condition.value.
    
    %% Properties
    properties
        value;
    end
    
    %% Methods
    methods
             
        function obj = Initial_Condition(value)
            
            obj.value = value; 
            
        end
        
    end
end

