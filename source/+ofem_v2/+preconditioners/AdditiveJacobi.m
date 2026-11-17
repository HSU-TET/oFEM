classdef AdditiveJacobi < handle
    %ADDITIVEJACOBI two-level additive preconditioner for H(curl) problems
    %
    %   B = sum_K R_K' (R_K A R_K')^-1 R_K  +  G L^-1 G'
    %
    %   The first part inverts the blocks of the assembled matrix A on the
    %   DOFs of every element (overlapping block Jacobi), the second part
    %   corrects on the gradient space with the discrete gradient G and
    %   L = G'*A*G.
    %
    %   Single matrix:
    %     p = AdditiveJacobi(dofs, A);
    %     p.nodalSolve(G, L, freeDOFs, amsDOFs);   % L = G'*A*G (all DOFs)
    %     y = p.preconditioner(r);
    %
    %   Frequency sweep for A = S + 1i*omega*M (eddy currents): the blocks of
    %   S and M are extracted once, and G'*S*G = 0, so L = 1i*omega*G'*M*G
    %   is factorised once for all frequencies:
    %     p = AdditiveJacobi(dofs, S, M);
    %     p.gradientCorrection(G, M, freeDOFs, amsDOFs, 'Solver', 'ict');
    %     for omega = ...
    %         p.setFrequency(omega);  % blocks of S + omega*M, L scaled by 1/(1i*omega)
    %         x = bicgstab(A, b, tol, maxit, @(r) p.preconditioner(r));
    %     end
    %   The block part uses S + omega*M (real, positive definite), the
    %   gradient part the exact L of A = S + 1i*omega*M.
    %
    %   gradientSolve(rhs) solves L*x = rhs accurately (e.g. to remove the
    %   gradient part of a right hand side).
    %
    %   preconditioner accepts blocks of right hand sides (one per column).

    properties
        prec;                   % sum of the block inverses on the free DOFs
        G;                      % discrete gradient (free DOFs x gradient DOFs)
        L;                      % factorisation of the gradient matrix
        nodalCorrect = 0;       % 1 if the gradient correction is used
        useDecomp = 0;          % 1: L is a decomposition, 0: L is an ICT factor
        gradScale = 1;          % L^-1 = gradScale * (factorised matrix)^-1
    end

    properties (Access = private)
        Ab;                     % blocks of A (single matrix), Ns x Ns x Ne
        Sb; Mb;                 % blocks of S and M (frequency sweep)
        rowF; colF; ic;         % sparsity pattern of the block sum on the free DOFs
        keep; nf;
        Lmat;                   % gradient matrix that is factorised
        sweep = false;
    end

    methods
        function obj = AdditiveJacobi(dofhandler,A,M)
            E = dofhandler.el2DOF;
            [Ne,Ns] = size(E);
            Nd = dofhandler.Nd;
            obj.sweep = nargin > 2;
            if obj.sweep
                C = A + 1i*M;           % extract S and M in one pass
            else
                C = A;
            end
            B = extractBlocks(C,E,Nd);
            if obj.sweep
                obj.Sb = real(B);
                obj.Mb = imag(B);
            else
                obj.Ab = B;
            end

            % sparsity pattern of the sum of the block inverses, restricted
            % to the free DOFs; reused for every frequency
            I = repmat(E,1,Ns)';
            J = repelem(E',Ns,1);
            [ukey,~,obj.ic] = unique((J(:)-1)*Nd + I(:));
            ri = mod(ukey-1,Nd)+1;
            ci = (ukey-ri)/Nd+1;
            freeMap = zeros(Nd,1);
            freeMap(dofhandler.freeDOFs) = 1:numel(dofhandler.freeDOFs);
            obj.keep = freeMap(ri) > 0 & freeMap(ci) > 0;
            obj.rowF = freeMap(ri(obj.keep));
            obj.colF = freeMap(ci(obj.keep));
            obj.nf = numel(dofhandler.freeDOFs);

            if ~obj.sweep
                obj.prec = obj.assembleInverse(obj.Ab);
                obj.Ab = [];
            end
        end

        function setFrequency(obj,omega)
            % block part for S + omega*M, gradient part for S + 1i*omega*M
            if ~obj.sweep
                error('ofem_v2:AdditiveJacobi:noSweep', ...
                    'setFrequency needs the constructor AdditiveJacobi(dofs,S,M).');
            end
            obj.prec = obj.assembleInverse(obj.Sb + omega*obj.Mb);
            obj.gradScale = 1/(1i*omega);
        end

        function nodalSolve(obj,G,L,freeDOFs,amsdofs)
            % gradient correction with a given L = G'*A*G (indexed by all
            % gradient DOFs)
            obj.nodalCorrect = 1;
            Gf = G(freeDOFs,amsdofs);
            nz = full(any(Gf,1));               % zero columns make L singular
            obj.G = Gf(:,nz);
            L = L(amsdofs,amsdofs);
            obj.Lmat = L(nz,nz);
            obj.gradScale = 1;
            obj.L = factorise(obj.Lmat);
            obj.useDecomp = 1;
        end

        function gradientCorrection(obj,G,M,freeDOFs,amsdofs,varargin)
            % gradient correction for A = S + 1i*omega*M with G'*S*G = 0:
            % L = 1i*omega*G'*M*G, factorised once ('chol' exact, 'ict'
            % incomplete Cholesky with 'DropTol', much cheaper to apply)
            ip = inputParser;
            ip.addParameter('Solver','chol',@(s) any(strcmpi(s,{'chol','ict'})));
            ip.addParameter('DropTol',1e-3);
            ip.parse(varargin{:});
            obj.nodalCorrect = 1;
            Gf = G(freeDOFs,amsdofs);
            obj.G = Gf(:,full(any(Gf,1)));      % zero columns make G'MG singular
            K = obj.G'*M(freeDOFs,freeDOFs)*obj.G;
            obj.Lmat = (K+K')/2;
            if strcmpi(ip.Results.Solver,'chol')
                obj.L = factorise(obj.Lmat);
                obj.useDecomp = 1;
            else
                obj.L = ictFactor(obj.Lmat,ip.Results.DropTol);
                obj.useDecomp = 0;
            end
        end

        function x = gradientSolve(obj,rhs,tol)
            % accurate solution of L*x = rhs
            if nargin < 3, tol = 1e-12; end
            if obj.useDecomp
                x = obj.L\rhs;
            else
                Lf = obj.L;
                x = zeros(size(rhs),'like',rhs);
                for k = 1:size(rhs,2)
                    [x(:,k),flag] = pcg(obj.Lmat,rhs(:,k),tol,1000,Lf,Lf');
                    if flag ~= 0
                        warning('ofem_v2:AdditiveJacobi:gradientSolve', ...
                            'pcg for the gradient system did not converge (flag %d).',flag);
                    end
                end
            end
            x = obj.gradScale*x;
        end

        function y = preconditioner(obj,r)
            if isempty(obj.prec)
                error('ofem_v2:AdditiveJacobi:noFrequency', ...
                    'Call setFrequency before using the preconditioner.');
            end
            y = obj.prec*r;
            if obj.nodalCorrect
                rn = obj.G'*r;
                if obj.useDecomp
                    e = obj.L\rn;
                else
                    e = obj.L'\(obj.L\rn);
                end
                y = y + obj.G*(obj.gradScale*e);
            end
        end
    end

    methods (Access = private)
        function P = assembleInverse(obj,B)
            vals = accumarray(obj.ic,reshape(pageinv(B),[],1));
            P = sparse(obj.rowF,obj.colF,vals(obj.keep),obj.nf,obj.nf);
        end
    end
end

function d = factorise(A)
    % Cholesky if Hermitian positive definite, LU otherwise
    d = [];
    if ishermitian(A)
        try
            d = decomposition(A,'chol','CheckCondition',false);
        catch
        end
    end
    if isempty(d)
        d = decomposition(A,'lu','CheckCondition',false);
    end
end

function L = ictFactor(A,droptol)
    alpha = 0;
    while true
        try
            L = ichol(A,struct('type','ict','droptol',droptol,'diagcomp',alpha));
            return
        catch err
            if ~strcmp(err.identifier,'MATLAB:ichol:Breakdown') || alpha > 1
                rethrow(err);
            end
            alpha = max(2*alpha,1e-4);
        end
    end
end

function B = extractBlocks(C,E,Nd)
    % B(:,:,e) = full(C(E(e,:),E(e,:))) for all elements at once. Indexing
    % the sparse matrix element by element costs O(size of C) per element;
    % here every entry is found by binary search in the sorted (column
    % major) keys of the nonzeros.
    [Ne,Ns] = size(E);
    [ri,ci,v] = find(C);
    skey = (ci-1)*Nd + ri;              % strictly increasing (CSC order)
    clear ri ci
    edges = [skey; Inf];
    B = zeros(Ns,Ns,Ne,'like',v);
    for j = 1:Ns
        for i = 1:Ns
            q = (E(:,j)-1)*Nd + E(:,i);
            loc = discretize(q,edges);
            hit = ~isnan(loc);
            hit(hit) = skey(loc(hit)) == q(hit);   % entry stored in C
            b = zeros(Ne,1,'like',v);
            b(hit) = v(loc(hit));
            B(i,j,:) = reshape(b,1,1,Ne);
        end
    end
end
