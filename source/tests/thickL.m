function [k] = thickL(deg)
%THICKL Summary of this function goes here
%   Detailed explanation goes here
    mesh = ofem_v2.Geometry;
    mesh.load_from_msh('thickL');
    mesh.create_edges();
    mesh.create_faces();

    mesh.connectFa2Ed;

    fe = ofem_v2.elements.NedelecZagl(mesh.dim,deg);
    fe.computeBasis();

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

    S = phys.S(dofs.freeDOFs,dofs.freeDOFs);
    M = phys.M(dofs.freeDOFs,dofs.freeDOFs);
    
    [R,P,C] = equilibrate(M);
    Sn = R*P*S*C;
    Mn = R*P*M*C;

    opt.useprecond = 5.5;

    opts.tol = 1e-14;

    %ksq = eigifp(S,M,200,opt);
    %ksq = irbleigs(S,M,opts);

    [v,ksq] = eigs(Sn,Mn,20,10,opts);
    %k2 = eigs(S,M,20,-1e7);
    v1 = zeros(length(dofs.DOFs),size(v,2));
    v1(dofs.freeDOFs,:) = v;
    %ksq = sort(sqrt(ksq));
    %ksq(abs(ksq)<1e-1) = [];
    %ksq = sort(ksq);
    ksq = diag(ksq);
    k = sort(ksq);
    k(k<1) = [];
    k = k(1:8);

    %reconstructTest
end

