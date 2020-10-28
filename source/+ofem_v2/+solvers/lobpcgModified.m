function [X,l,Wn,lh] = lobpcgModified(X,S,M,C,n,maxit,tol)
%MYLOBPCG Summary of this function goes here
%   Detailed explanation goes here
    
    tau = 10*eps;
    
    while norm(full(X'*M*X-eye(size(X,2)))) > tau
        X = svqb(M,X,10*eps);
    end

    m = size(X,2);
    [Cx,~,l] = rayleighRitz(X,S,eye(size(S)),m,0);
    X = X*Cx;
    
    while norm(full(X'*M*X-eye(size(X,2)))) > 10*eps
        X = svqb(M,X,10*eps);
    end
    
    J = 1:size(X,2);
    I = ones(size(J));
    
    Wn = zeros(maxit,size(X,2));
    lh = zeros(maxit,size(X,2));
    
    %[L,U] = ilu(S+sigma*M);
    L = ichol(S);
    q = 1;
    idx = logical(ones(size(X,2),1));
    useortho = 1;
    
    P = [];
    
    for it = 1:maxit
        % Compute residuals
        W = S*X-M*X*l;
        Wn(it,:) = vecnorm(W);
        lh(it,:) = diag(l);
        % New convergence
        Wn(it,:) = vecnorm(W)./((normest(S)+abs(diag(l))'.*normest(M)).*vecnorm(X));
        I(Wn(it,:)<tol) = 0;
        
        if size(X,2)-sum(I) >= n
            X = X(:,~logical(I));
            l = diag(l);
            l = l(~logical(I));
            disp('converged!')
            break;
        end
        % Soft locking in order
        if Wn(it,q) < tol
            idx(q) = 0;
            disp(['Eig ',num2str(q),' converged at iteration ',num2str(it)]);
            disp('removing from pool');
            q = q+1;
        end
        W = W(:,idx);

        for i = 1:size(W,2)
            [W(:,i),~] = pcg(S,W(:,i),1e-8,1000,L,L');
        end
        
        %W = (S)\W;
        W = orthoDrop(M,W,[X,P],tau);
                
        [Cx,Cp,l] = rayleighRitzImproved([X,P,W],S,M,m,sum(idx),useortho);
        
        %[~,i] = sort(abs(diag(l)));
        l = l(1:m,1:m);
        
        X = [X,P,W]*Cx;
        P = [X,P,W]*Cp;
        
        P = orthoDrop(M,P,X,tau);
        
        if it == maxit
%             X = X(:,~logical(I));
%             l = diag(l);
%             l = l(~logical(I));
            disp('maximum iterations reached!')
            disp([num2str(sum(I)),' eigenvalues did not converge'])
        end
        restart = 0;
    end
end

function U = ortho(M,U,V,tau)
    tau_r = 1e-10;
    while norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
        if norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
            disp(norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))))
        end
        U = U-V*(V'*M*U);
        while norm(full(U'*M*U-eye(size(U,2)))) > tau
            U = svqb(M,U,tau_r);
        end
    end
end

function U = orthoDrop(M,U,V,tau)
    tau_r = 10*eps;
    tau_drop = 10*eps;
    
    j = 0;
    
    while norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
        U = U-V*(V'*M*U);
        if norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
            disp(norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))))
        end
        while norm(full(U'*M*U-eye(size(U,2))))/(norm(full(M*U))*norm(full(U))) > tau
            j = j+1;
            if j==1
                [U,skip] = svqb(M,U,tau_r);
            else
                U = svqbDrop(M,U,tau_drop);
            end
        end
        if skip
            break;
        end
        j = 0;
    end
end

function [U,skip] = svqb(M,U,tau_r)
    D = diag(diag(U'*M*U))^(-1/2);
    [Z,theta] = svd(D*U'*M*U*D);
    theta = diag(theta);
    J = (theta < tau_r*max(max(theta)));
    theta(J) = tau_r*max(max(theta));
    theta = diag(theta);
    U = U*D*Z*theta^(-1/2);
    if cond(theta)<4
        skip = 1;
    else
        skip = 0;
    end
end

function U = svqbDrop(M,U,tau_drop)
    D = diag(diag(U'*M*U))^(-1/2);
    [Z,theta] = svd(D*U'*M*U*D);
    J = ~(diag(theta) < tau_drop*max(max(theta)));
    U = U*D*Z(:,J)*theta(J,J)^(-1/2);
end
    
function [Cx,Cp,theta] = rayleighRitz(C,S,M,m,q)
    D = diag(diag(C'*M*C))^(-1/2);
    R = chol(D*C'*M*C*D);
    [Z,theta] = svd(R'\D*C'*S*C*D/R);
    [~,i] = sort(diag(theta));
    theta = theta(i,i);
    Z = Z(:,i);
    Cx = D*R\Z(:,1:m);
    Cp = [];
    if q > 0
        [Q,~] = qr(Z(1:m,m+1:m+q)');
        Cp = D*R\Z(:,m+1:m+q)*Q';
    end
end

function [Cx,Cp,theta,useOrtho] = rayleighRitzImproved(C,S,M,m,q,useOrtho)
    if useOrtho
        [Z,theta] = svd(C'*S*C);
        [~,i] = sort(diag(theta));
        theta = theta(i,i);
        Z = Z(:,i);
        Cx = Z(:,1:m);
        %[Q,~] = qr(Z(1:m,m+1:m+q)');
        %Cp = Z(:,m+1:m+q)*Q';
        Cp = [zeros(size(Cx,2),size(Cx,2));Z(m+1:end,1:m)];
        
    else
        D = diag(diag(C'*M*C))^(-1/2);
        R = chol(D*C'*M*C*D);
        [Z,theta] = svd(R'\D*C'*S*C*D/R);
        [~,i] = sort(diag(theta));
        theta = theta(i,i);
        Z = Z(:,i);
        Cx = D*R\Z(:,1:m);
        [Q,~] = qr(Z(1:m,m+1:end)');
        Cp = D*R\Z(:,m+1:end)*Q';
        useOrtho = 1;
    end
    %theta = blkdiag(theta(1:m,1:m),Q(:,1:m-q)'*(theta(m+1:m+q,m+1:m+q))*Q(:,1:m-q));
        
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    