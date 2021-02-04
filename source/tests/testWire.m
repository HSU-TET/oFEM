close all
clear

mesh = ofem_v2.Geometry;
mesh.load_from_msh('testWire');

mesh.create_edges;
mesh.create_faces;

mesh.connectFa2Ed;

mesh.jacobiandata;

fe = ofem_v2.elements.NedelecZagl(3,2);
[N,curlN] = fe.computeBasis;

% fe = ofem_v2.elements.loadFE('NE2');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);
dofs.generateDOFs();

air = ofem_v2.materials.Material();
air.kappa = air.eps0;
air.mu = 1/air.mu0;

cond = ofem_v2.materials.Material;
cond.kappa = 5.6e7;
cond.mu = 1/cond.mu0;

PEC = ofem_v2.boundary.DirichletEdge([0;0;0],'PEC',mesh);

force = ofem_v2.Volume_force([0;0;1e6],'Conductor',mesh);

phys = ofem_v2.Physical_Problem(fe,mesh,1,1,0);
mesh.setMaterial('Conductor',cond);
mesh.setMaterial('Air',air);
phys.setParaS('mu');
phys.setParaM('kappa');
phys.setBoundaryCondition(PEC);
phys.attachDOFHandler(dofs);

phys.assemble();

S = phys.S;
M = phys.M;
b = phys.b;

free = dofs.freeDOFs;

f = 1e9;
w = 2*pi*f;

A = S+1i*w*M;

u = zeros(dofs.Nd,1);


%%

[L,U] = ilu(A(free,free));

%u(free) = bicgstab(A(free,free),b(free),1e-6,100,L,U);
u(free) = gmres(A(free,free),b(free),40,1e-12,100,L,U);
















