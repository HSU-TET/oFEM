%% Preconditioners for H(curl) systems under mesh refinement
%
% PCG for (S + M) u = b on the cube with PEC boundary (curl-curl plus mass,
% unit coefficients; the same type of operator as S + omega*M in the eddy
% current problem) with
%
%   none            plain CG
%   AdditiveJacobi  element-block Jacobi plus a correction on the gradients
%                   (exact Cholesky of G'*A*G)
%   HiptmairXu      auxiliary space preconditioner (hypre AMS type) with
%                   algebraic multigrid for the auxiliary H1 problems,
%                   multiplicative and additive cycle
%
% The iteration count of HiptmairXu stays constant when the mesh is
% refined, the others grow.

close all
clear

order = 0;              % H(curl) order
hs    = [2 4 6 8 10];   % mesh refinements
tol   = 1e-8;

fprintf('%3s %9s | %13s | %13s | %13s | %13s\n','h','free DOFs','none','AdditiveJacobi', ...
    'HX mult.','HX add.');
for h = hs
    %% Mesh and matrices
    file = sprintf('cube_h%d',h);
    if ~isfile([file,'.msh'])
        system(sprintf('gmsh cube.geo -setnumber h %d -3 -v 0 -o %s.msh',h,file));
    end
    mesh = ofem_v2.Geometry;
    mesh.load_from_msh(file);
    mesh.reorderAC;
    mesh.create_faces;
    mesh.create_edges;
    mesh.connectFa2Ed;
    mesh.jacobiandata;
    unit = ofem_v2.materials.Material;
    unit.stiff = 1;
    unit.mass = 1;
    mesh.setMaterial('Domain',unit);

    fe = ofem_v2.elements.loadFE(sprintf('HCurl_3D_Order_%d',order));
    dofs = ofem_v2.DOFHandler(mesh);
    dofs.attach(fe);
    phys = ofem_v2.Physical_Problem(fe,mesh,1,1,0);
    phys.setParaS('stiff');
    phys.setParaM('mass');
    phys.setBoundaryCondition(ofem_v2.boundary.DirichletEdge([0;0;0],'PEC',mesh));
    phys.attachDOFHandler(dofs);
    dofs.generateDOFs;
    phys.assemble(1);

    A = phys.S + phys.M;
    free = dofs.freeDOFs;
    Af = A(free,free);
    rng(1);
    b = Af*randn(numel(free),1);

    %% Preconditioners
    res = cell(1,4);
    % plain CG
    t = tic; [~,~,~,it] = pcg(Af,b,tol,5000); res{1} = [it,toc(t)];
    % AdditiveJacobi with exact gradient correction
    t = tic;
    aj = ofem_v2.preconditioners.AdditiveJacobi(dofs,A);
    G = fe.AMS(phys);
    aj.nodalSolve(G,G'*A*G,free,dofs.AMSDOFs(phys));
    [~,~,~,it] = pcg(Af,b,tol,5000,@(r) aj.preconditioner(r)); res{2} = [it,toc(t)];
    % Hiptmair-Xu, multiplicative and additive
    cycles = {'multiplicative','additive'};
    for c = 1:2
        t = tic;
        hx = ofem_v2.preconditioners.HiptmairXu.fromProblem(phys,A,'Cycle',cycles{c});
        [~,~,~,it] = pcg(Af,b,tol,5000,@(r) hx.apply(r)); res{2+c} = [it,toc(t)];
    end
    fprintf('%3d %9d |',h,numel(free));
    fprintf(' %4d it %5.2fs |',cell2mat(res));
    fprintf('\n');
end
