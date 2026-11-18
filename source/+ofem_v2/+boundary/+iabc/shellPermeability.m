function mur = shellPermeability(delta,order,dim)
%SHELLPERMEABILITY relative permeabilities of the improvised asymptotic
%boundary condition (IABC) shells, Meeker, IEEE Trans. Magn. 49(10), 2013
%
%   mur = ofem_v2.boundary.iabc.shellPermeability(delta, order)        3D / axisymmetric
%   mur = ofem_v2.boundary.iabc.shellPermeability(delta, order, 2)     2D planar
%
%   The domain of radius R is surrounded by ORDER shells of thickness
%   d = delta*R each, mur(k) is the relative permeability of shell k
%   (k = 1 innermost). The outer boundary of the last shell has B.n = 0
%   (n x A = 0, i.e. DirichletEdge for H(curl) elements). The shells
%   reproduce the impedance of free space exactly for the harmonics
%   n = 1..ORDER at r = R.
%
%   The matching conditions are solved with Newton's method, continued
%   in delta from Meeker's small delta limits (the system has several
%   roots; this selects the one of Tables I/II in the paper).
    if nargin < 3, dim = 3; end
    if delta <= 0, error('iabc:delta','delta = d/R has to be positive.'); end
    if order == 1 && dim == 3                 % closed form, eq. (16)
        mur = (delta^3+3*delta^2+3*delta+3)/(delta*(delta^2+3*delta+3));
        return
    elseif order == 1                         % eq. (10)
        mur = (delta^2+2*delta+2)/(delta*(delta+2));
        return
    end

    d0 = min(delta,1e-3);
    mur = smallDelta(d0,order,dim);
    steps = unique([d0, logspace(log10(d0),log10(delta),40), delta]);
    for s = steps
        mur = newton(@(m) residual(m,s,order,dim),mur);
    end
end

function mur = smallDelta(delta,order,dim)
    % small delta limits, first rows of Tables I and II
    switch 10*dim + order
        case 32, mur = [4*delta; 2/delta];
        case 33, mur = [1/(7*delta); 120*delta/49; 20/(7*delta)];
        case 22, mur = [3*delta; 3/(2*delta)];
        case 23, mur = [1/(6*delta); 10*delta/6; 10/(6*delta)];
        otherwise
            error('iabc:order','Orders 1-3 are supported.');
    end
end

function F = residual(mur,delta,order,dim)
    % admittance mismatch (relative to the free space value) for n = 1..order,
    % radius normalised to R = 1
    r = 1 + (0:order)*delta;
    F = zeros(order,1);
    for n = 1:order
        if dim == 3, Yfree = -(n+1); else, Yfree = -n; end
        F(n) = admittance(mur,r,n,dim)/Yfree - 1;
    end
end

function Y = admittance(mur,r,n,dim)
    % mu_r * dpsi/dr / psi at r(1) for the n-th harmonic of the magnetic
    % scalar potential, dpsi/dr = 0 (B.n = 0) at r(end)
    if dim == 3, p = -(n+1); else, p = -n; end      % psi = a r^n + b r^p
    Y = 0;
    for k = numel(mur):-1:1
        ro = r(k+1); ri = r(k); m = mur(k);
        ab = [ro^n, ro^p; m*n*ro^(n-1), m*p*ro^(p-1)] \ [1; Y];
        Y = m*(n*ab(1)*ri^(n-1) + p*ab(2)*ri^(p-1)) / (ab(1)*ri^n + ab(2)*ri^p);
    end
end

function x = newton(F,x)
    % Newton in log(mu) (keeps the permeabilities positive)
    y = log(x);
    G = @(y) F(exp(y));
    for it = 1:50
        g = G(y);
        if norm(g) < 1e-13, break, end
        J = zeros(numel(y));
        for k = 1:numel(y)
            e = zeros(size(y)); e(k) = 1e-7;
            J(:,k) = (G(y+e)-g)/1e-7;
        end
        y = y - J\g;
    end
    if norm(G(y)) > 1e-9
        error('iabc:newton','No convergence for the shell permeabilities.');
    end
    x = exp(y);
end
