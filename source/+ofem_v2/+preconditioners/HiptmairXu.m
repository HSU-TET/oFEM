classdef HiptmairXu < handle
    %HIPTMAIRXU auxiliary space preconditioner for H(curl) problems
    %(Hiptmair & Xu, SIAM J. Numer. Anal. 45, 2007; like hypre's AMS)
    %
    %   For A = curl-curl + mass (S + M, S + omega*M, ...) on H(curl):
    %
    %     B = smoother  +  G B_G G'  +  sum_d Pi_d B_d Pi_d'
    %
    %   G     discrete gradient (H1 -> H(curl)), B_G ~ (G'AG)^-1
    %   Pi_d  interpolation of the d-th component of continuous piecewise
    %         linear vector fields into the lowest-order (Whitney) edge DOFs,
    %         B_d ~ (Pi_d' A Pi_d)^-1
    %
    %   The auxiliary problems are scalar H1 problems, solved approximately
    %   with smoothed aggregation AMG (default), or with Cholesky / ICT.
    %   With AMG the number of iterations does not grow with the mesh.
    %   For elements of order >= 1 the higher-order modes are handled by the
    %   smoother (symmetric Gauss-Seidel on A).
    %
    %   hx = HiptmairXu(A, G, Pi)          A, G, Pi on the free DOFs,
    %                                      Pi = {Pi_x, Pi_y(, Pi_z)}
    %   hx = HiptmairXu.fromProblem(phys, A)   builds G and Pi (A: all DOFs)
    %   y = hx.apply(r)                    (also hx.preconditioner(r))
    %
    %   Options: 'Cycle'      'multiplicative' (default, symmetric:
    %                         smooth, G, Pi, G, smooth) or 'additive'
    %            'AuxSolver'  'amg' (default), 'chol', 'ict'
    %
    %   A has to be real symmetric positive definite. For complex problems
    %   like A = S + 1i*omega*M use the preconditioner of S + omega*M, it can
    %   be applied to complex residuals (e.g. with bicgstab or gmres).

    properties
        A;              % system matrix (free DOFs)
        G;              % discrete gradient (free DOFs x nodes)
        Pi;             % cell array of interpolation blocks
        cycle = 'multiplicative';
        auxSolver = 'amg';
    end

    properties (Access = private)
        BG;             % solver for G'AG
        BPi;            % solvers for Pi_d' A Pi_d
        L; U; Dg;       % Gauss-Seidel factors and diagonal of A
    end

    methods
        function obj = HiptmairXu(A,G,Pi,varargin)
            ip = inputParser;
            ip.addParameter('Cycle','multiplicative',@(s) any(strcmpi(s,{'multiplicative','additive'})));
            ip.addParameter('AuxSolver','amg',@(s) any(strcmpi(s,{'amg','chol','ict'})));
            ip.parse(varargin{:});
            obj.cycle = lower(ip.Results.Cycle);
            obj.auxSolver = lower(ip.Results.AuxSolver);
            if ~isreal(A)
                error('ofem_v2:HiptmairXu:complex', ...
                    'A has to be real (use e.g. S + omega*M for S + 1i*omega*M).');
            end
            if ~iscell(Pi), Pi = {Pi}; end

            obj.A = (A+A')/2;
            obj.G = dropZeroColumns(G);
            obj.Pi = cellfun(@dropZeroColumns,Pi,'UniformOutput',false);
            obj.L = tril(obj.A);
            obj.U = triu(obj.A);
            obj.Dg = full(diag(obj.A));

            obj.BG = obj.auxiliary(obj.G);
            obj.BPi = cellfun(@(P) obj.auxiliary(P),obj.Pi,'UniformOutput',false);
        end

        function x = apply(obj,r)
            A = obj.A;
            switch obj.cycle
                case 'multiplicative'
                    x = obj.L\r;                                    % forward GS
                    x = x + obj.correctG(r - A*x);
                    x = x + obj.correctPi(r - A*x);
                    x = x + obj.correctG(r - A*x);
                    x = x + obj.U\(r - A*x);                        % backward GS
                case 'additive'
                    % symmetric Gauss-Seidel + both corrections
                    x = obj.U\(obj.Dg.*(obj.L\r));
                    x = x + obj.correctG(r) + obj.correctPi(r);
            end
        end

        function x = preconditioner(obj,r)
            x = obj.apply(r);
        end
    end

    methods (Static)
        function hx = fromProblem(phys,A,varargin)
            % G and Pi from the H(curl) problem phys, A on all DOFs (e.g.
            % phys.S + phys.M); everything is restricted to the free DOFs
            [G,Pi] = ofem_v2.preconditioners.HiptmairXu.auxiliaryMatrices(phys);
            f = phys.DOFs.freeDOFs;
            % gradients of H1 functions that vanish on the Dirichlet
            % boundary (otherwise G'AG is singular)
            G = G(f,phys.DOFs.AMSDOFs(phys));
            Pi = cellfun(@(P) P(f,:),Pi,'UniformOutput',false);
            hx = ofem_v2.preconditioners.HiptmairXu(A(f,f),G,Pi,varargin{:});
        end

        function [G,Pi] = auxiliaryMatrices(phys)
            % G: discrete gradient of the H(curl) element (all DOFs x H1
            % DOFs, from HCurlElement.AMS); Pi{d}: interpolation of the
            % d-th component of piecewise linear vector fields at the mesh
            % nodes into the Whitney DOFs (edge e has DOF e).
            fe = phys.element;
            g = phys.geometry;
            G = fe.AMS(phys);
            dim = g.dim;
            Nd = phys.DOFs.Nd;
            ed = g.ed;
            Ne = size(ed,1);
            % Whitney DOF of a linear field u: circulation along the edge
            % ed(e,1) -> ed(e,2), t_e.(u_a + u_b)/2. The H(curl) basis uses
            % N_e = l_a grad l_b ... with the opposite sign of the textbook
            % Whitney function, hence the factor -1.
            t = reshape(g.co(:,1,ed(:,2)) - g.co(:,1,ed(:,1)),dim,[]);
            % components that are round-off (edges perpendicular to an
            % axis) would give numerically zero columns in Pi
            t(abs(t) < 1e-12*vecnorm(t,2,1)) = 0;
            rows = [1:Ne,1:Ne]';
            cols = [ed(:,1);ed(:,2)];
            Pi = cell(1,dim);
            for d = 1:dim
                vals = -[t(d,:),t(d,:)]'/2;
                Pi{d} = sparse(rows,cols,vals,Nd,g.Nco);
            end
        end
    end

    methods (Access = private)
        function B = auxiliary(obj,P)
            Ap = P'*obj.A*P;
            Ap = (Ap+Ap')/2;
            switch obj.auxSolver
                case 'amg'
                    amg = ofem_v2.preconditioners.SmoothedAggregationAMG(Ap);
                    B = @(r) amg.apply(r);
                case 'chol'
                    d = decomposition(Ap,'chol','CheckCondition',false);
                    B = @(r) d\r;
                case 'ict'
                    Lc = ictFactor(Ap);
                    Uc = Lc';
                    B = @(r) Uc\(Lc\r);
            end
        end

        function x = correctG(obj,r)
            x = obj.G*obj.BG(obj.G'*r);
        end

        function x = correctPi(obj,r)
            x = zeros(size(r),'like',r);
            for d = 1:numel(obj.Pi)
                x = x + obj.Pi{d}*obj.BPi{d}(obj.Pi{d}'*r);
            end
        end
    end
end

function P = dropZeroColumns(P)
    % drop (numerically) zero columns, they make the auxiliary matrix singular
    % (columns of Pi below 1e-8 come from round-off in nearly
    % perpendicular edges; genuine columns are > 1e-6 relative)
    n = full(sqrt(sum(P.^2,1)));
    P = P(:,n > 1e-8*max(n));
end

function L = ictFactor(A)
    alpha = 0;
    while true
        try
            L = ichol(A,struct('type','ict','droptol',1e-3,'diagcomp',alpha));
            return
        catch err
            if ~strcmp(err.identifier,'MATLAB:ichol:Breakdown') || alpha > 1
                rethrow(err);
            end
            alpha = max(2*alpha,1e-4);
        end
    end
end
