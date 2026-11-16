function [X,l,Wn,lh] = lobpcgModProj(X,S,M,~,n,maxit,tol,Y,YMY,dofs)
%LOBPCGMODPROJ LOBPCG for S*x = lambda*M*x with the gradient space
%projected out (Maxwell eigenvalue problems).
%
%   [X,l,Wn,lh] = lobpcgModProj(X,S,M,C,n,maxit,tol,Y,YMY,dofs)
%
%   X      initial block (all DOFs or only the free DOFs), m columns
%   S, M   stiffness and mass matrix (all DOFs)
%   C      unused, kept for compatibility. The preconditioner is (S+M)^-1.
%   n      number of wanted eigenpairs (n <= m)
%   maxit  maximum number of iterations
%   tol    tolerance for the relative backward error
%          ||S*x - lambda*M*x|| / ((||S|| + |lambda|*||M||)*||x||)
%   Y      discrete gradients, columns must vanish on the non-free DOFs
%   YMY    Y'*M*Y
%   dofs   free DOFs
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

    tau = 1e-14;

    Wn = zeros(maxit,m);
    lh = zeros(maxit,m);

    dYMY = spdFactor(YMY);
    dSM = spdFactor(S+M);
    project = @(U) U - Y*(dYMY\(MY*U));

    % Norm estimates for the backward error
    NS = svds(S,1,'largest','SubspaceDimension',60,'Tolerance',1e-10);
    NM = svds(M,1,'largest','SubspaceDimension',60,'Tolerance',1e-10);

    % Initial block: project out gradients, M-orthonormalise, Rayleigh-Ritz
    X = orthonormalize(M,X,10*eps,tau);
    X = project(X);
    X = orthonormalize(M,X,10*eps,tau);
    [Z,l] = ritz(X'*S*X);
    X = X*Z;

    P = [];
    q = 1;                 % first not yet converged eigenpair
    active = true(1,m);    % columns whose residual is still expanded
    converged = false;

    for it = 1:maxit
        % Residuals and backward errors
        W = S*X - (M*X).*l';
        lh(it,:) = l;
        Wn(it,:) = vecnorm(W)./((NS+abs(l').*NM).*vecnorm(X));

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
            P = P(:,active);
        end

        X = orthonormalize(M,X,tau,tau);
        if ~isempty(P)
            P = orthoDrop(M,P,X,tau);
        end

        % Precondition, project out gradients, orthogonalise against [X,P]
        W = project(dSM\W);
        W = orthoDrop(M,W,[X,P],tau);

        % Rayleigh-Ritz on the M-orthonormal basis [X,W,P]
        B = [X,W,P];
        [Z,theta] = ritz(B'*S*B);
        Cx = Z(:,1:m);
        l = theta(1:m);
        X = B*Cx;
        P = B(:,m+1:end)*Cx(m+1:end,:);
    end

    if ~converged
        warning('lobpcgModProj:maxit', ...
            'Maximum iterations reached, %d of %d eigenpairs did not converge.', ...
            n-q+1,n);
    end
end

function d = spdFactor(A)
    % Cholesky factorisation of a symmetric positive definite matrix,
    % LU if it is not positive definite
    A = (A+A')/2;
    try
        d = decomposition(A,'chol','CheckCondition',false);
    catch
        d = decomposition(A,'lu','CheckCondition',false);
    end
end

function [Z,theta] = ritz(G)
    % Eigenpairs of the symmetric projected matrix, ascending
    [Z,T] = eig((G+G')/2);
    [theta,i] = sort(diag(T));
    Z = Z(:,i);
end

% The orthogonalisation routines carry MU = M*U along, so M is applied
% once per Gram-Schmidt pass, and compute 2-norms of tall matrices from
% their small Gram matrices.

function n = norm2(A)
    n = sqrt(norm(A'*A));
end

function e = orthoError(U,MU)
    e = norm(U'*MU-eye(size(U,2)))/(norm2(MU)*norm2(U));
end

function U = orthonormalize(M,U,tau_r,tau)
    % Repeated SVQB until U'*M*U = I
    MU = M*U;
    for k = 1:10
        if orthoError(U,MU) <= tau
            return
        end
        [U,MU] = svqb(U,MU,tau_r,false);
    end
    warning('lobpcgModProj:ortho','M-orthonormalisation did not reach %g.',tau);
end

function U = orthoDrop(M,U,V,tau)
    % M-orthogonalise U against V and M-orthonormalise U, dropping
    % linearly dependent columns
    tau_drop = 10*eps;
    MV = M*V;
    nMV = norm2(MV);
    for k = 1:3
        U = U-V*(MV'*U);
        MU = M*U;
        for j = 1:10
            if isempty(U) || orthoError(U,MU) <= tau
                break
            end
            [U,MU] = svqb(U,MU,tau_drop,j > 1);
        end
        if isempty(U) || norm(MV'*U)/(nMV*norm2(U)) <= tau
            return
        end
    end
end

function [U,MU] = svqb(U,MU,tau,drop)
    % SVQB step. Small eigenvalues of the Gram matrix are either clamped
    % to tau*max (drop = false) or their directions dropped (drop = true).
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
    U = U*T;
    MU = MU*T;
end
