classdef  Material
    % Material class specifies the material properties 
    
    properties  %(SetAccess= private)
       lambda; 
       rho;         % density; 
       epsilon;     % permittivity;
       mu;          % permeability;
       e_mod;       % e_modulus;
       kappa;       % thermal_conductivity;
       c;           % specific_heat_capacity;
	   stiff;
	   mass;
        
	end
	
	properties(Constant)
		eps0 = 8.8541878176E-12;
		mu0 = 4*pi*1e-7;
	end
    
    methods
         function obj = material(lambda, rho, epsilon, my, e_mod, kappa, c)
            % Material class specifies the material properties 
            obj.lambda = lambda; 
            obj.rho = rho;
            obj.epsilon = epsilon; 
            obj.my = my; 
            obj.e_mod = e_mod;
            obj.kappa = kappa; 
            obj.c = c; 
         end
         
         function setMaterial(obj, physicalProblem, partName)
           for i = 1:size(physicalProblem.geometry.parts,1)  
             if isequal(physicalProblem.geometry.parts{1,i} ,partName)
                 physicalProblem.geometry.parts{4,i} = obj; 
                 break; 
             end 
           end
         end
          

        
         
    end
end

