function [X,l,Wn,lh] = LOBPCGBasic(X,S,M,C,maxit,tol,Y,dofs)
%MYLOBPCG Summary of this function goes here
%   Detailed explanation goes here

    % M-Orthonormalize X
    XMX = X'*M*X;
    if ~isreal(XMX)
        XMX = (XMX+XMX')/2;
    end
    R = chol(XMX);
    X = X/R;
        
    condestHist = zeros(maxit,1);
    condestHist(1) = 1;
    
    explicit = 0;
    
    % Compute initial Ritz Values
    XSX = X'*S*X;
    XSX = (XSX+XSX')/2;
    [v,l] = eig(XSX);
    X = X*v;
    
    J = 1:size(X,2);
    I = ones(size(J));
    
    Wn = zeros(maxit,size(X,2));
    lh = zeros(maxit,size(X,2));
    
    setup.droptol = 1e-9;
    %setup.milu = 'row';
    setup.udiag = 1;
    setup.type = 'nofill';
    
    sigma = 1.9;
    restart = 0;
    
    %[L,U] = ilu(S+sigma*M);
    L = ichol(S+2*M);
    q = 1;
    idx = logical(ones(size(X,2),1));
    
    m = size(X,2);
    
    for it = 1:maxit
        % Compute residuals
        W = S*X-M*X*l;
        Wn(it,:) = vecnorm(W);
        lh(it,:) = diag(l);
%         explicit = 1;
        if any(Wn(it,:)<1e-4)||it == 1
            explicit = 1;
        else
            explicit = 0;
        end
        I(Wn(it,:)<tol) = 0;
        
        if sum(I) <= 1
            X = X(:,~logical(I));
            l = diag(l);
            l = l(~logical(I));
            disp('converged!')
            break;
        end
        % Soft locking in order
        if Wn(it,q) < tol
            idx(q) = 0;
            q = q+1;
            disp(['Eig ',num2str(q),' converged at iteration ',num2str(it)]);
            disp('removing from pool');
        end
        W = W(:,idx);

        % Apply preconditioner
        for i = 1:size(W,2)
            [W(:,i),~] = pcg(S+2*M,W(:,i),1e-8,1000,L,L');
        end
%         W = (S)\W;
        
        % Orthonormalize W to X
        W = W-X*((M*X)'*W);

        % M-Orthonormalize W
        WMW = W'*M*W;
        WMW = (WMW+WMW')/2;
        [R,f] = chol(WMW);
        if f==0
            W = W/R;
        else
            warning('myLOBPCG:ResidualNotFullRank','Residual does not have full rank!')
            break;
        end
        
        condestMean = mean(max(1,it-10-round(log(sum(idx)))):it);
        
        for j = 1:3
            if it > 1 && f == 0
                % M-Orthonormalize P
                P = P;
                PMP = P'*M*P;
                PMP = (PMP+PMP')/2;
                [R,f] = chol(PMP);
                if f==0 && ~restart
                    P = P/R;
                    if ~ explicit
%                         Compute Gram matrices
                        gramS = [l        , X'*S*W,  zeros(m);
                                (X'*S*W)', W'*S*W,  W'*S*P;
                                zeros(m),(W'*S*P)',P'*S*P];
                        gramM = [eye(m),  X'*M*W,         zeros(m,sum(idx));
                                 (X'*M*W)',       eye(sum(idx)), W'*M*P;
                                 zeros(sum(idx),m),       (W'*M*P)',      eye(sum(idx))];
                    else
                        % explicit gram
                        gramS = [X'*S*X  , X'*S*W,  X'*S*P;
                                (X'*S*W)', W'*S*W,  W'*S*P;
                                (X'*S*P)',(W'*S*P)',P'*S*P];
                        gramM = [X'*M*X   , X'*M*W,  X'*M*P;
                                (X'*M*W)',  W'*M*W,  W'*M*P;
                                (X'*M*P)', (W'*M*P)',P'*M*P];
                    end
                else
                    if ~explicit
                        gramS = [l,        X'*S*W;
                                (X'*S*W)',W'*S*W];
                        gramM = [eye(size(X,2)), X'*M*W;
                                (X'*M*W)',      eye(size(X,2))];
                    else
                    % explicit gram
                        gramS = [X'*S*X,  X'*S*W;
                                (X'*S*W)',W'*S*W];
                        gramM = [X'*M*X,  X'*M*W;
                                (X'*M*W)',W'*M*W];
                    end
                end


            else
                % Gram matrices for initial run
                if ~explicit
                    gramS = [l,        X'*S*W;
                             (X'*S*W)',W'*S*W];
                    gramM = [eye(size(X,2)), X'*M*W;
                             (X'*M*W)',      eye(size(X,2))];
                else
                    % explicit gram
                    gramS = [X'*S*X,  X'*S*W;
                            (X'*S*W)',W'*S*W];
                    gramM = [X'*M*X,  X'*M*W;
                            (X'*M*W)',W'*M*W];
                end
            end
            
            condestHist(it) = log10(condest(gramM))+1;
            if (condestHist(it)/condestMean > 2 && condestHist(it) > 10 )|| condestHist(it) > 8
                if j == 1 && ~restart
                    restart=1; %steepest descent restart for stability
                elseif j==3
                    warning('myLOBPCG:Condition',...
                    'Gramm matrix ill-conditioned: results unpredictable');
                    break;
                else
                    warning('myLOBPCG:Condition',...
                    'Gramm matrix ill-conditioned: setting explicit');
                    explicit = 1;
                end
            else
                break;
            end
            
            
            
        end
        
        % Compute Ritz values
        [V,l] = eig(gramS,gramM);
        
        if it > 1 && f == 0 && ~restart
            [~,i] = sort(abs(diag(l)));
            l = l(i(1:m),i(1:m));
            Vx = V(1:m,i(1:m));
            Vw = V(m+1:m+sum(idx),i(1:m));
            Vp = V(m+sum(idx)+1:m+2*sum(idx),i(1:m));
            
            P = W*Vw+P*Vp;
            X = X*Vx+P;
            
        else
            [l,i] = sort(abs(diag(l)));
            l = diag(l(1:size(X,2)));
            Vx = V(1:m,i(1:m));
            Vw = V(m+1:m+sum(idx),i(1:m));
            P = W*Vw;
            X = X*Vx+P;

        end
        
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
    while normest(V'*M*U)/(normest(M*U)*normest(U)) < tau
        U = U-V*(V'*M*U);
        while normest(U'*M*U-eye(size(U,2))) < tau
            U = svqb(M,U,tau_r);
        end
    end
end

function U = orthoDrop(M,U,V,tau)
    tau_r = 1e-10;
    tau_drop = 1e-10;
    
    j = 0;
    
    while normest(V'*M*U)/(normest(M*U)*normest(U)) < tau
        U = U-V*(V'*M*U);
        while normest(U'*M*U-eye(size(U,2))) < tau
            j = j+1;
            if j==1
                U = svqb(M,U,tau_r);
            else
                U = svqbDrop(M,U,tau_drop);
            end
        end
        j = 0;
    end
end

function U = svqb(M,U,tau_r)
    D = diag(U'*M*U)^(-1/2);
    [Z,theta] = svd(D*U'*M*U*D);
    t = [];
    for i = diag(theta)
        if i < tau_r*max(max(theta))
            t = [t,tau_r*max(max(theta))]
        else
            t = [t,i];
        end
    end
    U = U*D*Z*diag(t)^(-1/2);
end

function U = svqbDrop(M,U,tau_drop)
    D = diag(U'*M*U)^(-1/2);
    [Z,theta] = svd(D*U'*M*U*D);
    J = ones(size(Z,2),1);
    J(diag(theta) < tau_drop*max(max(theta))) = 0;
    U = U*D*Z(:,J)*theta(J,J)^(-1/2);
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    