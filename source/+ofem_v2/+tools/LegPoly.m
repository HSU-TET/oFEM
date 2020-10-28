classdef LegPoly < handle
	%LEGENDRE Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		leg = []
		legS = []
		legInt = []
		legIntS = []
        legIntNB = []
	end
	
	methods
		function obj = LegPoly()
			%obj.order = order;
		end
		
		function computePolynomials(obj,order)
			syms x t c
			obj.leg = [1,x];
			obj.legS = [1,x];
			obj.legInt = [x,(x^2-1)/2];
			obj.legIntS = [x,(x^2-t^2)/2];
            obj.legIntNB = [-0.5*c,-0.5*x*c];
			for n = 1:order
				obj.leg(n+2) = ((2*n+3)*obj.leg(n+1)*x-(n+1)*obj.leg(n))/(n+2);
				obj.legS(n+2) = ((2*n+3)*x*obj.legS(n+1)-(n+1)*t^2*obj.legS(n))/(n+2);
            end
            for n = 2:order
				obj.legInt(n+1) = ((2*n-1)*x*obj.legInt(n)-(n-2)*obj.legInt(n-1))/(n+1);
				obj.legIntS(n+1) = ((2*n-1)*x*obj.legIntS(n)-(n-2)*t^2*obj.legIntS(n-1))/(n+1);
            end
            for n = 2:order
                obj.legIntNB(n+1) = (2*(n+2)-3)/n*obj.legIntNB(n)*x-(n-3)/n*t^2*obj.legIntNB(n-1);
            end
		end
	end
end

