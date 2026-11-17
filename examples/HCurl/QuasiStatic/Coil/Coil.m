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
c
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
Lgpu = gpuArray(L);
Ugpu = gpuArray(U);
Agpu = gpuArray(complex(A(dofsNe.freeDOFs,dofsNe.freeDOFs)));
bgpu = gpuArray(b(dofsNe.freeDOFs));


%%
tic
A_vec = v_potential.u;
A_vecgpu = gpuArray(A_vec(dofsNe.freeDOFs));
Agpu = gpuArray(complex(A(dofsNe.freeDOFs,dofsNe.freeDOFs)));
bgpu = gpuArray(b(dofsNe.freeDOFs));
toc
tic
A_vecgpu = bicgstab(Agpu,bgpu,1e-6,2000);
toc
tic
A_vecgpu = bicgstab(Agpu,bgpu,1e-6,5000);
toc
A_vec(dofsNe.freeDOFs) = gather(A_vecgpu);
toc

%%
tic
[~,A_rec] = ofem_v2.tools.reconstruct(v_potential,full(A_vec),1);
toc
[~,B] = ofem_v2.tools.reconstructCurl(v_potential,full(A_vec));
toc

%%

mesh.export_UCD([pwd,'/export'], 'field', {'B_i',imag(B)','','Cell'},{'A_i',imag(A_rec)','','Cell'},...
										  {'B_r',real(B)','','Cell'},{'A_r',real(A_rec)','','Cell'},...
										  {'B_a',abs(B)','','Cell'},{'A_a',abs(A_rec)','','Cell'});









