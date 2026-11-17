classdef SmoothedAggregationAMG < handle
    %SMOOTHEDAGGREGATIONAMG algebraic multigrid for symmetric positive
    %definite (Laplace/mass type) matrices
    %
    %   amg = SmoothedAggregationAMG(A)
    %   amg = SmoothedAggregationAMG(A,'Theta',0.08,'MaxCoarse',2000)
    %   x = amg.apply(b)          % one V-cycle, b may have several columns
    %
    %   Smoothed aggregation (Vanek, Mandel, Brezina): strength of
    %   connection |a_ij| >= theta*sqrt(a_ii*a_jj), greedy aggregation,
    %   piecewise constant tentative prolongation smoothed with one damped
    %   Jacobi step, Galerkin coarse matrices. The V-cycle uses a forward
    %   Gauss-Seidel sweep before and a backward sweep after the coarse
    %   correction, so it is a symmetric preconditioner (for PCG). The
    %   coarsest level is solved with a Cholesky factorisation.

    properties
        levels;                 % struct array: A, P, L, U
        coarse;                 % factorisation on the coarsest level
        theta = 0.08;
        maxCoarse = 2000;
        maxLevels = 15;
    end

    methods
        function obj = SmoothedAggregationAMG(A,varargin)
            ip = inputParser;
            ip.addParameter('Theta',obj.theta);
            ip.addParameter('MaxCoarse',obj.maxCoarse);
            ip.addParameter('MaxLevels',obj.maxLevels);
            ip.parse(varargin{:});
            obj.theta = ip.Results.Theta;
            obj.maxCoarse = ip.Results.MaxCoarse;
            obj.maxLevels = ip.Results.MaxLevels;

            A = (A+A')/2;
            lv = struct('A',{},'P',{},'L',{},'U',{});
            for k = 1:obj.maxLevels
                n = size(A,1);
                lv(k).A = A;
                lv(k).L = tril(A);
                lv(k).U = triu(A);
                if n <= obj.maxCoarse || k == obj.maxLevels
                    break
                end
                agg = aggregate(A,obj.theta);
                nc = max(agg);
                if nc >= n || nc == 0
                    break                       % no coarsening possible
                end
                % tentative prolongation: normalised indicator of aggregates
                cnt = accumarray(agg,1);
                T = sparse((1:n)',agg,1./sqrt(cnt(agg)),n,nc);
                % smoothing with damped Jacobi, omega = 4/(3 rho(D^-1 A))
                d = full(diag(A));
                omega = 4/(3*spectralRadius(A,d));
                P = T - omega*(spdiags(1./d,0,n,n)*(A*T));
                lv(k).P = P;
                A = P'*A*P;
                A = (A+A')/2;
            end
            obj.levels = lv;
            Ac = lv(end).A;
            try
                obj.coarse = decomposition(Ac,'chol','CheckCondition',false);
            catch
                obj.coarse = decomposition(Ac,'lu','CheckCondition',false);
            end
        end

        function x = apply(obj,b)
            x = obj.vcycle(1,b);
        end

        function c = complexity(obj)
            % operator complexity sum(nnz(A_k))/nnz(A_1)
            c = sum(arrayfun(@(l) nnz(l.A),obj.levels))/nnz(obj.levels(1).A);
        end
    end

    methods (Access = private)
        function x = vcycle(obj,k,b)
            lv = obj.levels(k);
            if k == numel(obj.levels)
                x = obj.coarse\b;
                return
            end
            x = lv.L\b;                                 % forward Gauss-Seidel
            x = x + lv.P*obj.vcycle(k+1,lv.P'*(b-lv.A*x));
            x = x + lv.U\(b-lv.A*x);                    % backward Gauss-Seidel
        end
    end
end

function agg = aggregate(A,theta)
    % greedy aggregation on the strength graph
    n = size(A,1);
    d = full(abs(diag(A)));
    [i,j,v] = find(A);
    strong = i ~= j & abs(v) >= theta*sqrt(d(i).*d(j));
    S = sparse(i(strong),j(strong),true,n,n);
    S = S | S';
    [ri,cj] = find(S);                                  % column-wise neighbour lists
    ptr = [0; cumsum(accumarray(cj,1,[n,1]))];
    nb = @(c) ri(ptr(c)+1:ptr(c+1));

    agg = zeros(n,1);
    k = 0;
    % phase 1: nodes whose neighbourhood is still free start an aggregate
    for c = 1:n
        if agg(c) == 0
            N = nb(c);
            if all(agg(N) == 0)
                k = k+1;
                agg([c;N]) = k;
            end
        end
    end
    % phase 2: attach the remaining nodes to a neighbouring aggregate
    rest = find(agg == 0);
    agg1 = agg;
    for c = rest'
        N = nb(c);
        a = agg1(N);
        a = a(a > 0);
        if ~isempty(a)
            agg(c) = a(1);
        end
    end
    % phase 3: new aggregates from what is left (incl. isolated nodes)
    for c = find(agg == 0)'
        if agg(c) == 0
            N = nb(c);
            N = N(agg(N) == 0);
            k = k+1;
            agg([c;N]) = k;
        end
    end
end

function rho = spectralRadius(A,d)
    % power iteration for rho(D^-1 A)
    x = rand(size(A,1),1);
    rho = 1;
    for it = 1:15
        y = (A*x)./d;
        rho = norm(y)/norm(x);
        x = y/norm(y);
    end
end
