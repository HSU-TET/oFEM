function res = movingConductor(meshFile,varargin)
%MOVINGCONDUCTOR magnetic field of a conductor moving through a uniform
%background field, open boundary emulated with improvised asymptotic
%boundary conditions (IABC, Meeker 2013)
%
%   res = movingConductor(meshFile, 'Name', value, ...)
%
%   Solves  curl(1/mu curl A) + eps*A - kappa*v x curl A = 0  (lowest-order
%   H(curl) elements) with the background field B0 imposed by the IABC
%   shells (ofem_v2.boundary.iabc.uniformField) and n x A = 0 on the outer boundary.
%
%   Mesh: physical volumes Conductor (name below), 'Air', 'Layer1' ...
%   'LayerN' (concentric spherical shells of equal thickness), physical
%   surface 'Outer' (outer boundary of the last shell).
%
%   Options
%     'Conductor'       part name of the conductor ('Coil'), '' for none
%     'Mur'             relative permeability of the conductor (1)
%     'Kappa'           conductivity of the conductor (0)
%     'Velocity'        velocity of the conductor (zeros(3,1))
%     'B0'              background field (zeros(3,1))
%     'Source'          'motion'     velocity term kappa*v x curl A (default)
%                       'firstOrder' no velocity term, source kappa*v x B0
%                                    in the conductor (no background field)
%     'Regularization'  mass term eps of all materials (1)
%     'ShellMur'        IABC permeabilities ([] computed from the mesh)
%     'Solver'          'gmres' (default) or 'qmr'
%     'Preconditioner'  'hx' (default) Hiptmair-Xu of S + M (ofem_v2) as a
%                       right preconditioner, or 'jacobi'
%     'Solve'           false: return after the assembly, res.system holds
%                       the system on the free DOFs (A, b)
%     'Tol', 'MaxIt'    solver tolerance (1e-5) and iterations (2000). The
%                       attainable residual is limited by rounding in the
%                       gradient part of A, which scales with
%                       kappa*|v|/Regularization (copper, kappa*|v| = 5.8e7:
%                       2e-5 for eps = 1, 2e-7 for eps = 100 with B changed
%                       by 1e-5); mesh refinement does not lower it
%     'Restart'         gmres restart (100)
%
%   res: mesh, phys, dofs, u, X (element centres), A, B, Bind = B - B0,
%   J = kappa*v x B, Jind = kappa*v x Bind (in the conductor), resvec
%   (relative residuals of the solver), residual (true relative residual), flag, mur (IABC shells), R, d, times
    ip = inputParser;
    ip.addParameter('Conductor','Coil');
    ip.addParameter('Mur',1);
    ip.addParameter('Kappa',0);
    ip.addParameter('Velocity',zeros(3,1));
    ip.addParameter('B0',zeros(3,1));
    ip.addParameter('Source','motion',@(s) any(strcmpi(s,{'motion','firstOrder'})));
    ip.addParameter('Regularization',1);
    ip.addParameter('ShellMur',[]);
    ip.addParameter('Solver','gmres',@(s) any(strcmpi(s,{'qmr','gmres'})));
    ip.addParameter('Tol',1e-5);
    ip.addParameter('MaxIt',2000);
    ip.addParameter('Restart',100);
    ip.addParameter('Preconditioner','hx',@(s) any(strcmpi(s,{'hx','jacobi'})));
    ip.addParameter('Solve',true);
    ip.parse(varargin{:});
    o = ip.Results;
    v = o.Velocity(:);
    B0 = o.B0(:);
    firstOrder = strcmpi(o.Source,'firstOrder');
    tStart = tic;

    [~,~,ext] = fileparts(meshFile);
    if isempty(ext), meshFile = [meshFile '.msh']; end
    if ~isfile(meshFile)
        geo = [meshFile(1:end-4) '.geo'];
        if ~isfile(geo)
            error('movingConductor:mesh','Neither %s nor %s found.',meshFile,geo);
        end
        fprintf('creating %s with gmsh\n',meshFile);
        [st,out] = system(sprintf('gmsh "%s" -3 -v 0 -o "%s"',geo,meshFile));
        if st ~= 0, error('movingConductor:gmsh','gmsh failed:\n%s',out); end
    end
    mesh = ofem_v2.Geometry;
    mesh.load_from_msh(meshFile);
    fe = ofem_v2.elements.loadFE('HCurl_3D_Order_0');
    mu0 = 4e-7*pi;

    % materials
    air = ofem_v2.materials.Material;
    air.stiff = 1/mu0;
    air.mass = o.Regularization;
    mesh.setMaterial('Air',air);
    hasConductor = ~isempty(o.Conductor);
    if hasConductor
        cond = ofem_v2.materials.Material;
        cond.stiff = 1/(o.Mur*mu0);
        cond.mass = o.Regularization;
        cond.kappa = o.Kappa;
        mesh.setMaterial(o.Conductor,cond);
        cIdx = mesh.parts{3,strcmp(mesh.parts(1,:),o.Conductor)};
    end

    % open boundary and background field
    layers = ofem_v2.boundary.iabc.layerNames(mesh);
    [mur,R,d] = ofem_v2.boundary.iabc.addShells(mesh,'Layers',layers,'Mass',o.Regularization,'Mur',o.ShellMur);
    if firstOrder
        mesh.setForce(o.Conductor,ofem_v2.Volume_force(o.Kappa*cross(v,B0),o.Conductor,mesh));
    elseif any(B0)
        ofem_v2.boundary.iabc.uniformField(mesh,B0,layers);
    end

    % assembly
    dofs = ofem_v2.DOFHandler(mesh);
    dofs.attach(fe);
    phys = ofem_v2.Physical_Problem(fe,mesh,1,1,0);
    phys.setParaS('stiff');
    phys.setParaM('mass');
    phys.setBoundaryCondition(ofem_v2.boundary.DirichletEdge([0;0;0],'Outer',mesh));
    phys.attachDOFHandler(dofs);
    dofs.generateDOFs;
    phys.assemble(1);
    A = phys.S + phys.M;
    if hasConductor && ~firstOrder && any(v) && o.Kappa ~= 0
        A = A + fe.assembleVelocityTerm(phys,cIdx,o.Kappa,v);
    end
    res.times.assembly = toc(tStart);

    free = dofs.freeDOFs;
    Af = A(free,free);
    b = phys.b(free);
    if ~o.Solve
        res.system = struct('A',Af,'b',b,'free',free);
        res.mesh = mesh; res.phys = phys; res.dofs = dofs;
        res.mur = mur; res.R = R; res.d = d;
        return
    end

    % preconditioner: Hiptmair-Xu of the symmetric positive definite part
    % S + M, or Jacobi. The velocity term D makes A nonsymmetric, and it
    % only acts on the test side (D*G = 0, but G'*D is large compared to
    % G'*M*G), so HX is applied as a right preconditioner,
    % A*P^-1*y = b, x = P^-1*y: the Krylov method then sees the true
    % residual. As a left preconditioner HX amplifies G'*D*x by orders of
    % magnitude and gmres stagnates.
    t = tic;
    switch lower(o.Preconditioner)
        case 'hx'
            hx = ofem_v2.preconditioners.HiptmairXu.fromProblem(phys,phys.S + phys.M);
            Pinv = @(r) hx.apply(r);
            op = @(y,tr) rightOperator(Af,Pinv,y,tr);   % hx is symmetric
            P = [];
        case 'jacobi'
            P = spdiags(full(diag(Af)),0,numel(free),numel(free));
            op = Af;
    end
    res.times.preconditioner = toc(t);

    % solution
    t = tic;
    switch lower(o.Solver)
        case 'qmr'
            if isempty(P)
                [x,res.flag,~,~,resvec] = qmr(op,b,o.Tol,o.MaxIt);
            else
                [x,res.flag,~,~,resvec] = qmr(op,b,o.Tol,o.MaxIt,P);
            end
        case 'gmres'
            if isempty(P)
                [x,res.flag,~,~,resvec] = gmres(@(y) op(y,'notransp'),b,o.Restart,o.Tol, ...
                    ceil(o.MaxIt/o.Restart));
            else
                [x,res.flag,~,~,resvec] = gmres(op,b,o.Restart,o.Tol, ...
                    ceil(o.MaxIt/o.Restart),P);
            end
    end
    if isempty(P), x = Pinv(x); end
    res.times.solve = toc(t);
    if res.flag ~= 0
        warning('movingConductor:solver','%s did not converge (flag %d).',o.Solver,res.flag);
    end
    u = phys.u;
    u(free) = x;
    res.resvec = resvec/norm(b);
    res.residual = norm(b - Af*x)/norm(b);    % true relative residual

    % fields per element
    [X,Arec] = ofem_v2.tools.reconstruct(phys,full(u),1);
    [~,B] = ofem_v2.tools.reconstructCurl(phys,full(u));
    B = double(B);
    if firstOrder, Bbg = zeros(3,1); else, Bbg = B0; end
    Bind = B - Bbg;
    J = zeros(size(B));
    Jind = J;
    if hasConductor
        J(:,cIdx) = o.Kappa*cross(repmat(v,1,numel(cIdx)),B(:,cIdx));
        Jind(:,cIdx) = o.Kappa*cross(repmat(v,1,numel(cIdx)),Bind(:,cIdx));
    end

    res.mesh = mesh; res.phys = phys; res.dofs = dofs; res.u = u;
    res.X = double(X); res.A = double(Arec); res.B = B; res.Bind = Bind;
    res.J = J; res.Jind = Jind;
    res.mur = mur; res.R = R; res.d = d;
    res.times.total = toc(tStart);
end

function z = rightOperator(A,Pinv,y,tr)
    % A*P^-1*y, and (A*P^-1)'*y = P^-1*A'*y for the symmetric P
    if strcmp(tr,'transp')
        z = Pinv(A'*y);
    else
        z = A*Pinv(y);
    end
end
