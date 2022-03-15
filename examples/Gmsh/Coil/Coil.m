close all;
clear;

file = 'Coil';

mesh = ofem_v2.Geometry;
mesh.load_from_msh(file);
%mesh.reorderAC;
mesh.co = mesh.co*1e-3;
mesh.create_faces;
mesh.create_edges;
mesh.connectFa2Ed;

%%
fe = ofem_v2.elements.loadFE('H1_3D_Order_1');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);

copper = ofem_v2.materials.Material;
copper.kappa = 5.96e7; % Wikipedia
copper.mu = copper.mu0; % Wikipedia

ferrite = ofem_v2.materials.Material;
ferrite.kappa = 1e7; % Wikipedia: Iron
ferrite.mu = 5000 * ferrite.mu0; % Wikipedia: Iron

air = ofem_v2.materials.Material;
air.kappa = 1e-12; % Small conductivity to stabilize the computation realistic would be 1e-12 to 1e-15
air.mu = 1;

mesh.setMaterial('Air', air);
mesh.setMaterial('Coil', copper);
mesh.setMaterial('Ferrite', ferrite);

%%
upperVoltage = 1;
lowerVoltage = 0;

voltage = ofem_v2.Physical_Problem(fe,mesh,1,0,0);
voltage.setParaS('kappa');
voltage.attachDOFHandler(dofs);

upperBD = ofem_v2.boundary.Dirichlet(upperVoltage,'Upper',mesh);
lowerBD = ofem_v2.boundary.Dirichlet(lowerVoltage,'Lower',mesh);

voltage.setBoundaryCondition(upperBD);
voltage.setBoundaryCondition(lowerBD);

dofs.generateDOFs;

voltage.assemble;
%voltage.solve;

%%
L = ichol(voltage.S(dofs.freeDOFs,dofs.freeDOFs));

u = voltage.u;
u(dofs.freeDOFs) = pcg(voltage.S(dofs.freeDOFs,dofs.freeDOFs),voltage.b(dofs.freeDOFs),1e-9,1000,L,L');
%%
E = -voltage.gradCell(full(u));

kappa = ofem_v2.tools.matrixarray(zeros(1,1,mesh.Nint));
kappa(1,1,mesh.parts{3,1}) = ferrite.kappa;
kappa(1,1,mesh.parts{3,2}) = copper.kappa;

J = kappa*E;

mesh.export_UCD([pwd,'/export'],'voltage',{'U',u,''},{'E',squeeze(E)','','Cell'}, {'J',squeeze(J)','','Cell'});

%%
feNe = ofem_v2.elements.loadFE('HCurl_3D_Order_0');

dofsNe = ofem_v2.DOFHandler(mesh);
dofsNe.attach(feNe);

v_potential = ofem_v2.Physical_Problem(feNe,mesh,1,1,0);
v_potential.attachDOFHandler(dofsNe);

v_potential.setParaS('kappa');
v_potential.setParaM('mu');

upperBD = ofem_v2.boundary.DirichletEdge([0,0,0]','Upper',mesh);
lowerBD = ofem_v2.boundary.DirichletEdge([0,0,0]','Lower',mesh);

v_potential.setBoundaryCondition(upperBD);
v_potential.setBoundaryCondition(lowerBD);

force = ofem_v2.Volume_force(J,'Coil',mesh);

mesh.setForce('Coil',force);

dofsNe.generateDOFs;
v_potential.assemble();

%%
S = v_potential.S;
M = v_potential.M;

b = v_potential.b;

omega = 2*pi*1000; % 1kHz

A = S + 1i*omega*M;

%L = ichol(M(dofsNe.freeDOFs,dofsNe.freeDOFs));
[L,U] = ilu(A(dofsNe.freeDOFs,dofsNe.freeDOFs));

A = bicgstab(A(dofsNe.freeDOFs,dofsNe.freeDOFs),b(dofsNe.freeDOFs),1e-6,1000,L,U);













