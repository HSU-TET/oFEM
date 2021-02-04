function [x,err] = pcgTest(A,b,prec,dofs,varargin)
%PCGTEST Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 5
        it = varargin{1};
        tol = 1e-6;
    elseif nargin == 6
        it = varargin{1};
        tol = varargin{2};
    else
        it = 100;
        tol = 1e-6;
    end
    
    x = zeros(size(b));
    
    b2n = norm(b(dofs));
    
    r = b;
    err = zeros(it+1,1);
    err(1) = b2n;
    
    for i = 1:it
        z = prec(A,r);
        if i == 1
            p = z(dofs);
            rho = z(dofs).'*r(dofs);
        else
            rho1 = rho;
            rho = z(dofs).'*r(dofs);
            beta = (rho)/(rho1);
            p = z(dofs)+beta*p;
        end
        q = A(dofs,dofs)*p;
        alpha = (rho)/(p.'*q);
        x(dofs) = x(dofs)+alpha*p;
        r(dofs) = r(dofs)-alpha*q;
        err(i+1) = norm(r(dofs));
        disp(['iteration ',num2str(i),' tolerance ',num2str(err(i)/b2n),'\n'])
        if err(i)/b2n < tol
            break;
        end
    end
end

