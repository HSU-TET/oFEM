function u = eqsTimeStepping(Sk,Se,free,uD,u,dt,n,callback)
%EQSTIMESTEPPING implicit Euler for the electro-quasistatic problem
%
%   div((kappa + epsilon d/dt) grad phi) = 0
%
%   u = eqsTimeStepping(Sk,Se,free,uD,u,dt,n,callback) performs n steps
%
%       (Sk + Se/dt) u_i = Se/dt u_{i-1},   u_i = uD(t_i) on the fixed DOFs
%
%   Sk, Se    stiffness matrices with conductivity and permittivity
%   free      free DOFs
%   uD        function handle, uD(t) returns a full vector with the
%             Dirichlet values on the fixed DOFs (other entries ignored)
%   u         initial potential (full vector)
%   callback  optional function handle callback(i,t,u) called after
%             every step, e.g. for the export
    if nargin < 8, callback = []; end
    A = Sk + Se/dt;
    Aff = A(free,free);
    try
        dA = decomposition((Aff+Aff')/2,'chol');
    catch
        dA = decomposition(Aff,'lu');
    end
    isFree = false(size(u));
    isFree(free) = true;
    for i = 1:n
        t = i*dt;
        uOld = u;
        u = uD(t);
        u(isFree) = 0;
        rhs = Se*uOld/dt - A*u;
        u(free) = dA\rhs(free);
        if ~isempty(callback)
            callback(i,t,u);
        end
    end
end
