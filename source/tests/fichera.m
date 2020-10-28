%% Compute eigenmodes of -1x1\0x1 3D fichera domain

function [k] = fichera(deg)
    mesh = ofem_v2.Geometry;
    mesh.load_from_msh('fichera');
    mesh.reorderAC();
    mesh.create_edges();
    mesh.create_faces();

    mesh.connectFa2Ed;

    fe = ofem_v2.elements.loadFE(['NE',num2str(deg)]);
%     fe = ofem_v2.elements.NedelecZagl(mesh.dim,deg);
%     fe.computeBasis();

    dofs = ofem_v2.DOFHandler(mesh);
    dofs.attach(fe);
    dofs.generateDOFs();

    mat = ofem_v2.materials.Material();
    mat.kappa = 1;

    PEC = ofem_v2.boundary.DirichletEdge([0;0;0],'Boundary',mesh);

    phys = ofem_v2.Physical_Problem(fe,mesh,1,1,0);
    mesh.setMaterial('Domain',mat);
    phys.setParaS('kappa');
    phys.setParaM('kappa');
    phys.setBoundaryCondition(PEC);
    phys.attachDOFHandler(dofs);

    phys.assemble();
    
    X0 = randn(length(dofs.freeDOFs),10);

    [S_g,Y_g] = fe.AMS(phys);


    S = phys.S;%(dofs.freeDOFs,dofs.freeDOFs);
    M = phys.M;%(dofs.freeDOFs,dofs.freeDOFs);

    wDofs = dofs.AMSDOFs(phys);

    tic
    [V,l,Wn,lh] = ofem_v2.solvers.lobpcgModProj(X0,S,M,speye(size(S)),3,300,1e-8,Y_g(:,wDofs),S_g(wDofs,wDofs),dofs.freeDOFs);
    toc
    k = l(1:3)';
%     ksq = diag(l);
%     rounded = round(ksq);
%     k = sort(ksq);

end