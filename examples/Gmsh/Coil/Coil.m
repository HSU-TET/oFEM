close all;
clear;

file = 'Coil';

mesh = ofem_v2.Geometry;
mesh.load_from_msh(file);
mesh.reorderAC;
mesh.co = mesh.co*1e-3;
mesh.jacobiandata;
mesh.create_faces;
mesh.create_edges;
mesh.connectFa2Ed;

%%
fe = ofem_v2.elements.loadFE('H1_3D_Order_1');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);

copper = ofem_v2.materials.Material;
copper.kappa = 5.96e7; % Wikipedia
copper.mu = 1/copper.mu0; % Wikipedia

ferrite = ofem_v2.materials.Material;
ferrite.kappa = 1e7; % Wikipedia: Iron
ferrite.mu = 1/(5000 * ferrite.mu0); % Wikipedia: Iron

air = ofem_v2.materials.Material;
air.kappa = 1; % Small conductivity to stabilize the computation realistic would be 1e-12 to 1e-15
air.mu = 1/air.mu0;

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
u(dofs.freeDOFs) = pcg(voltage.S(dofs.freeDOFs,dofs.freeDOFs),voltage.b(dofs.freeDOFs),1e-9,10000,L,L');
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

v_potential.setParaS('mu');
v_potential.setParaM('kappa');

upperBD = ofem_v2.boundary.DirichletEdge([0,0,0]','Upper',mesh);
lowerBD = ofem_v2.boundary.DirichletEdge([0,0,0]','Lower',mesh);
PECBD = ofem_v2.boundary.DirichletEdge([0,0,0]','PEC',mesh);

v_potential.setBoundaryCondition(upperBD);
v_potential.setBoundaryCondition(lowerBD);
v_potential.setBoundaryCondition(PECBD);

force = ofem_v2.Volume_force(J,'Coil',mesh);

mesh.setForce('Coil',force);

dofsNe.generateDOFs;
v_potential.assemble();

%%
S = v_potential.S;
M = v_potential.M;

b = v_potential.b;

omega = 2*pi*100000; % 1kHz

A = S + 1i*omega*M;
opts = [];
opts.type = 'crout';
opts.droptol = 1e-4;
%opts.michol = 'on';
%opts.diagcomp = 1e4;


%L = ichol(A(dofsNe.freeDOFs,dofsNe.freeDOFs),opts);
[L,U] = ilu(A(dofsNe.freeDOFs,dofsNe.freeDOFs));
disp('here!')

%%
A_vec = v_potential.u;
A_vec(dofsNe.freeDOFs) = bicgstab(A(dofsNe.freeDOFs,dofsNe.freeDOFs),b(dofsNe.freeDOFs),1e-6,20000,L,U);

%%

[~,A_rec] = ofem_v2.tools.reconstruct(v_potential,full(A_vec),1);
[~,B] = ofem_v2.tools.reconstructCurl(v_potential,full(A_vec));

%%

mesh.export_UCD([pwd,'/export'], 'field', {'B_i',imag(B)','','Cell'},{'A_i',imag(A_rec)','','Cell'},...
										  {'B_r',real(B)','','Cell'},{'A_r',real(A_rec)','','Cell'},...
										  {'B_a',abs(B)','','Cell'},{'A_a',abs(A_rec)','','Cell'});









