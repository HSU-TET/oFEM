function [X,l,Wn,lh] = lobpcgModProj(X,S,M,~,n,maxit,tol,Y,YMY,dofs,varargin)
%LOBPCGMODPROJ LOBPCG for S*x = lambda*M*x with the gradient space
%projected out (Maxwell eigenvalue problems).
%
%   [X,l,Wn,lh] = lobpcgModProj(X,S,M,C,n,maxit,tol,Y,YMY,dofs)
%   [X,l,Wn,lh] = lobpcgModProj(...,'Preconditioner','ict','DropTol',1e-3)
%
%   X      initial block (all DOFs or only the free DOFs), m columns
%   S, M   stiffness and mass matrix (all DOFs)
%   C      unused, kept for compatibility
%   n      number of wanted eigenpairs (n <= m)
%   maxit  maximum number of iterations
%   tol    tolerance for the relative backward error
%          ||S*x - lambda*M*x|| / ((||S|| + |lambda|*||M||)*||x||)
%   Y      discrete gradients, columns must vanish on the non-free DOFs.
%          Zero columns are removed.
%   YMY    Y'*M*Y
%   dofs   free DOFs
%
%   Options (name-value pairs)
%   'Preconditioner'  'chol' (default): exact Cholesky solve with S+M.
%                     'ict': one application of an incomplete Cholesky
%                     factor of S+M. Needs far less memory than 'chol' but
%                     more iterations; meant for problems where the
%                     Cholesky factor of S+M does not fit into memory.
%   'DropTol'         drop tolerance of the incomplete Cholesky factor
%                     (default 1e-3)
%   The gradient projection always uses a Cholesky factor of Y'*M*Y.
%
%   X      M-orthonormal Ritz vectors on the free DOFs (m columns)
%   l      Ritz values as column vector (m entries)
%   Wn     backward error per iteration and vector (maxit x m)
%   lh     Ritz values per iteration (maxit x m)
%
%   Converged eigenpairs are soft locked in order: their residuals are
%   no longer expanded, but the vectors stay in the Rayleigh-Ritz step.
%   The search directions are M-orthogonalised against the gradients
%   x <- x - Y*(Y'*M*Y)^-1*Y'*M*x.
%   The products S*X, M*X, S*P and M*P are carried through the iteration
%   and recomputed every 10 iterations and before convergence is accepted.

    ip = inputParser;
    ip.addParameter('Preconditioner','chol', ...
        @(s) any(strcmpi(s,{'chol','ict'})));
    ip.addParameter('DropTol',1e-3,@(x) isscalar(x) && x > 0);
    ip.parse(varargin{:});
    opts = ip.Results;

    if size(X,1) ~= numel(dofs)
        X = X(dofs,:);
    end
    m = size(X,2);
    if n > m
        error('lobpcgModProj:blockSize', ...
            'Number of wanted eigenpairs (%d) exceeds the block size (%d).',n,m);
    end

    MY = Y'*M(:,dofs);
    M = M(dofs,dofs);
    S = S(dofs,dofs);
    Y = Y(dofs,:);

    % Remove zero gradient columns, they make Y'*M*Y singular
    keep = full(any(Y,1));
    Y = Y(:,keep);
    MY = MY(keep,:);
    YMY = YMY(keep,keep);

    tau = 1e-14;
    refresh = 10;

    Wn = zeros(maxit,m);
    lh = zeros(maxit,m);

    dYMY = spdFactor(YMY);
    project = @(U) U - Y*(dYMY\(MY*U));
    switch lower(opts.Preconditioner)
        case 'chol'
            dSM = spdFactor(S+M);
            precond = @(R) dSM\R;
        case 'ict'
            L = ictFactor(sym2(S+M),opts.DropTol);
            U = L';
            precond = @(R) U\(L\R);
    end

    % Norm estimates for the backward error (only scale the tolerance)
    NS = svds(S,1,'largest','Tolerance',1e-4);
    NM = svds(M,1,'largest','Tolerance',1e-4);

    % Initial block: project out gradients, M-orthonormalise, Rayleigh-Ritz
    X = orthonormalize(M,X,10*eps,tau);
    X = project(X);
    X = orthonormalize(M,X,10*eps,tau);
    SX = S*X;
    MX = M*X;
    [Z,l] = ritz(X'*SX);
    X = X*Z; SX = SX*Z; MX = MX*Z;

    P = []; SP = []; MP = [];
    q = 1;                 % first not yet converged eigenpair
    active = true(1,m);    % columns whose residual is still expanded
    converged = false;
    lastRefresh = 0;

    for it = 1:maxit
        % Residuals and backward errors
        W = SX - MX.*l';
        Wn(it,:) = vecnorm(W)./((NS+abs(l').*NM).*vecnorm(X));

        % Accept convergence only with freshly computed products
        if q <= n && Wn(it,q) < tol && lastRefresh < it
            SX = S*X; MX = M*X; lastRefresh = it;
            W = SX - MX.*l';
            Wn(it,:) = vecnorm(W)./((NS+abs(l').*NM).*vecnorm(X));
        end
        lh(it,:) = l;

        % Soft locking in order
        while q <= n && Wn(it,q) < tol
            active(q) = false;
            fprintf('Eig %d converged at iteration %d\n',q,it);
            q = q+1;
        end
        if q > n
            converged = true;
            fprintf('All %d eigenpairs converged after %d iterations.\n',n,it);
            break;
        end

        W = W(:,active);
        if ~isempty(P)
            P = P(:,active); SP = SP(:,active); MP = MP(:,active);
        end

        % Refresh the carried products
        if it - lastRefresh >= refresh
            SX = S*X; MX = M*X;
            if ~isempty(P)
                SP = S*P; MP = M*P;
            end
            lastRefresh = it;
        end

        % X stays M-orthonormal; reorthonormalise only if it drifted
        [X,SX,MX] = orthonormalizeCarried(X,SX,MX,tau,tau);
        if ~isempty(P)
            [P,SP,MP] = orthoDropCarried(P,SP,MP,X,SX,MX,tau);
        end

        % Precondition, project out gradients, orthogonalise against [X,P]
        W = project(precond(W));
        [W,MW] = orthoDrop(M,W,[X,P],[MX,MP],tau);
        SW = S*W;

        % Rayleigh-Ritz on the M-orthonormal basis [X,W,P]
        B = [X,W,P]; SB = [SX,SW,SP]; MB = [MX,MW,MP];
        [Z,theta] = ritz(B'*SB);
        Cx = Z(:,1:m);
        l = theta(1:m);
        X = B*Cx; SX = SB*Cx; MX = MB*Cx;
        Cp = Cx(m+1:end,:);
        P = B(:,m+1:end)*Cp; SP = SB(:,m+1:end)*Cp; MP = MB(:,m+1:end)*Cp;
    end

    if ~converged
        warning('lobpcgModProj:maxit', ...
            'Maximum iterations reached, %d of %d eigenpairs did not converge.', ...
            n-q+1,n);
    end
end

function A = sym2(A)
    A = (A+A')/2;
end

function d = spdFactor(A)
    % Cholesky factorisation of a symmetric positive definite matrix,
    % LU if it is not positive definite
    A = sym2(A);
    try
        d = decomposition(A,'chol','CheckCondition',false);
    catch
        warning('lobpcgModProj:notSPD', ...
            'Matrix is not positive definite, using LU instead of Cholesky.');
        d = decomposition(A,'lu','CheckCondition',false);
    end
end

function L = ictFactor(A,droptol)
    % Threshold incomplete Cholesky, with diagonal compensation if a
    % nonpositive pivot occurs
    alpha = 0;
    while true
        try
            L = ichol(A,struct('type','ict','droptol',droptol,'diagcomp',alpha));
            return
        catch err
            if ~strcmp(err.identifier,'MATLAB:ichol:Breakdown') || alpha > 1
                rethrow(err);
            end
            alpha = max(2*alpha,1e-3);
        end
    end
end

function [Z,theta] = ritz(G)
    % Eigenpairs of the symmetric projected matrix, ascending
    [Z,T] = eig((G+G')/2);
    [theta,i] = sort(diag(T));
    Z = Z(:,i);
end

% 2-norms of tall matrices are computed from their small Gram matrices.

function n = norm2(A)
    n = sqrt(norm(A'*A));
end

function e = orthoError(U,MU)
    e = norm(U'*MU-eye(size(U,2)))/(norm2(MU)*norm2(U));
end

function U = orthonormalize(M,U,tau_r,tau)
    % M-orthonormalise U
    U = svqbLoop(U,M*U,[],tau_r,tau,false);
end

function [U,SU,MU] = orthonormalizeCarried(U,SU,MU,tau_r,tau)
    % M-orthonormalise U with carried products S*U and M*U
    [U,MU,SU] = svqbLoop(U,MU,SU,tau_r,tau,false);
end

function [U,MU] = orthoDrop(M,U,V,MV,tau)
    % M-orthogonalise U against V and M-orthonormalise U, dropping
    % linearly dependent columns. M*U is recomputed after every
    % projection because U may cancel strongly.
    tau_drop = 10*eps;
    nMV = norm2(MV);
    MU = [];
    for k = 1:3
        U = U-V*(MV'*U);
        MU = M*U;
        [U,MU] = svqbLoop(U,MU,[],tau_drop,tau,true);
        if isempty(U) || norm(MV'*U)/(nMV*norm2(U)) <= tau
            return
        end
    end
end

function [U,SU,MU] = orthoDropCarried(U,SU,MU,V,SV,MV,tau)
    % orthoDrop with carried products of U and V
    tau_drop = 10*eps;
    nMV = norm2(MV);
    for k = 1:3
        C = MV'*U;
        U = U-V*C; SU = SU-SV*C; MU = MU-MV*C;
        [U,MU,SU] = svqbLoop(U,MU,SU,tau_drop,tau,true);
        if isempty(U) || norm(MV'*U)/(nMV*norm2(U)) <= tau
            return
        end
    end
end

function [U,MU,SU] = svqbLoop(U,MU,SU,tau_r,tau,allowDrop)
    % Repeated SVQB until U'*M*U = I up to tau, or until a pass no
    % longer halves the error (round-off floor). The first pass clamps
    % small Gram eigenvalues, later passes drop them if allowDrop.
    e = orthoError(U,MU);
    for j = 1:10
        if isempty(U) || e <= tau
            return
        end
        T = svqb(U,MU,tau_r,allowDrop && j > 1);
        U = U*T; MU = MU*T;
        if ~isempty(SU), SU = SU*T; end
        eOld = e;
        e = orthoError(U,MU);
        if e > eOld/2
            break
        end
    end
    if e > 1e-10
        warning('lobpcgModProj:ortho','M-orthonormalisation only reached %.1e.',e);
    end
end

function T = svqb(U,MU,tau,drop)
    % SVQB transformation U -> U*T. Small eigenvalues of the Gram matrix
    % are either clamped to tau*max (drop = false) or their directions
    % dropped (drop = true).
    G = U'*MU;
    D = diag(diag(G).^(-1/2));
    G = D*G*D;
    [Z,theta] = eig((G+G')/2,'vector');
    if drop
        J = theta >= tau*max(abs(theta));
        Z = Z(:,J);
        theta = theta(J);
    else
        theta = max(theta,tau*max(abs(theta)));
    end
    T = D*Z*diag(theta.^(-1/2));
end
