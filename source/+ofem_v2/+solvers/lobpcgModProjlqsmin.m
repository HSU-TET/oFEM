function [X,l,Wn,lh] = lobpcgModProj(X,S,M,C,n,maxit,tol,Y,YMY,dofs)
%MYLOBPCG Summary of this function goes here
%   Detailed explanation goes here
    
	X = X(dofs,:);

    MY = (Y'*M);
    M = M(dofs,dofs);
    S = S(dofs,dofs);
    Y = Y(dofs,:);
    MY = MY(:,dofs);

	testTol1 = 1e-6;
	testTol2 = 1e-4;

    %L = ilu(YMY,struct('udiag',1,'droptol',1e-4));
    
    tau = 1e-14;
    
    J = 1:size(X,2);
    I = ones(size(J));
    
    Wn = zeros(maxit,size(X,2));
    lh = zeros(maxit,size(X,2));
    ltest = zeros(maxit,size(X,2));
    
    %[L,U] = ilu(S+sigma*M);
%     L1 = ichol(S+M);
    %L2 = ichol(YMY);
    %dYMY = decomposition(YMY,'lu','CheckCondition',false);
    %dSM = decomposition(S+M,'lu','CheckCondition',false);
	dYMY = pinv(YMY,1e-10);
	dSM = pinv(S+M,1e-10);
    L2 = [];
    q = 1;
    idx = logical(ones(size(X,2),1));
    useortho = 1;
    
    P = [];
    
    while norm(full(X'*M*X-speye(size(X,2))))/(norm(full(M*X))*norm(full(X))) > tau
        X = svqb(M,X,10*eps);
    end
    
    for j=1:1
%         disp(vecnorm(MY'*X))
        SX = zeros(size(YMY,1),size(X,2));
        %W = W-X(:,1:q-1)*(X(:,1:q-1)'*M*W);
%         for i = 1:size(X,2)
%             %temp = zeros(totalD,1);
%             temp = X(:,i);
%             [SX(:,i),~] = bicgstab(YMY,MY*temp,1e-6,1000,L2,L2');
%         end
        %SX = dYMY\(MY*X);
		SX = lsqminnorm(YMY,MY*X,testTol1);
        SX = Y*SX;
%         SX = SX(dofs,:);
%         SW = YMY\((MY')*W);
        X = X - SX;
    end
    
    NS = svds(S,1,'largest','SubspaceDimension',60,'Tolerance',1e-10);
    NM = svds(M,1,'largest','SubspaceDimension',60,'Tolerance',1e-10);

    while norm(full(X'*M*X-speye(size(X,2))))/(norm(full(M*X))*norm(full(X))) > tau
        X = svqb(M,X,10*eps);
    end


    m = size(X,2);
    [Cx,~,l] = rayleighRitz(X,S,speye(size(S)),m,0);
    X = X*Cx;
    
    SX = zeros(size(YMY,1),size(X,2));
    
%     for i = 1:size(X,2)
%         %temp = zeros(totalD,1);
%         temp = X(:,i);
%         [SX(:,i),~] = bicgstab(YMY,MY*temp,1e-6,1000,L2,L2');
%     end
    %SX = dYMY\(MY*X);
	%SX = lsqminnorm(YMY,MY*X,testTol1);
	SX = dYMY*(MY*X);
    SX = Y*SX;
    X = X - SX;

    while norm(full(X'*M*X-speye(size(X,2))))/(norm(full(M*X))*norm(full(X))) > tau
        X = svqb(M,X,tau);
    end
    
    for it = 1:maxit
        % Compute residuals
        
        W = S*X-M*X*l;
        ltest(it,:) = diag(X'*S*X);
        lh(it,:) = diag(l);
        % New convergence
        Wn(it,:) = vecnorm(W)./((NS+abs(diag(l))'.*NM).*vecnorm(X));
        I(Wn(it,:)<tol) = 0;
        
        

        % Soft locking in order
        while Wn(it,q) < tol
            idx(q) = 0;
            disp(['Eig ',num2str(q),' converged at iteration ',num2str(it)]);
            disp('removing from pool');
            P(:,1) = [];
            q = q+1;
            if q==n+1
                break;
            end
        end
        
        if q >= n+1
            I = zeros(size(J));
            X = X(:,~logical(I));
            l = diag(l);
            l = l(~logical(I));
            disp('converged!')
            break;
        end
        W = W(:,idx);
        
        while norm(full(X'*M*X-speye(size(X,2))))/(norm(full(M*X))*norm(full(X))) > tau
            X = svqb(M,X,tau);
        end
        if it>1
            P = orthoDrop(M,P,X,tau);
        end
        
%         for i = 1:size(W,2)
%             [W(:,i),~] = pcg(S+1.8*M,W(:,i),1e-8,1000,L1,L1');
%         end
		W = lsqminnorm(S+M,W,testTol1);
        %W = dSM\W;
        %W = (S+M)\W;
%         
        for j =1:1
            SW = zeros(size(YMY,1),size(W,2));
%             parfor i = 1:size(W,2)
%                 %temp = zeros(totalD,1);
%                 %[SW(:,i)] = gmres(YMY,MY*W(:,i),10,1e-4,1);
%                 SW(:,i) = decomp\(MY*W(:,i));
%             end
            %SW = dYMY\(MY*W);
			SW = lsqminnorm(YMY,MY*W,testTol2);
            SW = Y*SW;
%             SW = SW(dofs,:);
    %         SW = YMY\((MY')*W);
            W = W - SW;
        end
        
        W = orthoDrop(M,W,[X,P],tau);
        
        orthoOld = useortho;
        [Cx,Cp,l,useortho] = rayleighRitzImproved([X,W,P],S,M,m,sum(idx),useortho);
        
        if useortho == orthoOld
            l = l(1:m,1:m);

            X = [X,W,P]*Cx;
            %P = [X,W,P]*Cp;
            P = [W,P]*Cx(m+1:end,:);
            
        else
            l = diag(lh(it,:));
            continue;
        end
        
        if it == maxit
            disp('maximum iterations reached!')
            disp([num2str(sum(I)),' eigenvalues did not converge'])
        end
    end
end

function U = ortho(M,U,V,tau)
    tau_r = 1e-10;
    while norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
        if norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
            disp(norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))))
        end
        U = U-V*(V'*M*U);
        while norm(full(U'*M*U-speye(size(U,2)))) > tau
            U = svqb(M,U,tau_r);
        end
    end
end

function U = orthoDrop(M,U,V,tau)
    tau_r = 10*eps;
    tau_drop = 10*eps;
    
    j = 0;
    skip = 0;

    n1 = V'*M;
    n2 = norm(full(M*V));
    
    %while norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
    while norm(full(n1*U))/(n2*norm(full(U))) > tau
		U = U-V*(V'*M*U);
%         if norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))) > tau
%             disp(norm(full(V'*M*U))/(norm(full(M*V))*norm(full(U))))
%         end
        while norm(full(U'*M*U-speye(size(U,2))))/(norm(full(M*U))*norm(full(U))) > tau
            j = j+1;
            if j==1 %|| j==2
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
        [~,i] = sort(abs(diag(theta)));
        theta = theta(i,i);
        Z = Z(:,i);
        Cx = Z(:,1:m);
        %[Q,~] = qr(Z(1:m,m+1:m+q)');
        %Cp = Z(:,m+1:m+q)*Q';
        Cp = Z(m+1:end,1:m);
        
    else
        D = diag(diag(C'*M*C))^(-1/2);
        [R,f] = chol(D*C'*M*C*D);
        if f==0
            %[Z,theta] = svd(R'\D*C'*S*C*D/R);
            [Z,theta] = eig(D*C'*S*C*D,D*C'*M*C*D);
            [~,i] = sort(diag(theta));
            theta = theta(i,i);
            Z = Z(:,i);
            Cx = (D*R)\Z(:,1:m);
            [Q,L] = qr(Z(1:m,m+1:end)');
            Cp = (D*R)\Z(:,m+1:end)*Q(1:m,:)';
%             Cp = (D*R)\Z(:,m+1:m+q);
%             Cp = Cp(m+1:end,:);
        else
            useOrtho = 1;
            Cx = [];
            Cp = [];
            theta = [];
        end
    end
    %theta = blkdiag(theta(1:m,1:m),Q(:,1:m-q)'*(theta(m+1:m+q,m+1:m+q))*Q(:,1:m-q));
        
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    