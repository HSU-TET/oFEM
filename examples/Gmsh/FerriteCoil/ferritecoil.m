close all;
clear;

file = 'ferritecoil';
tic;
disp("Importing and setting up mesh.")
mesh = ofem_v2.Geometry;
mesh.load_from_msh(file);
mesh.reorderAC;
mesh.co = mesh.co*1e-3;
mesh.jacobiandata;
mesh.create_faces;
mesh.create_edges;
mesh.connectFa2Ed;
t = toc;
disp("Done in: ")
disp(t)

%%
tic;
disp("Setting up and assembling H1 Problem.")
fe = ofem_v2.elements.loadFE('H1_3D_Order_1');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);

copper = ofem_v2.materials.Material;
copper.kappa = 5.96e7; % Wikipedia
copper.mu = 1/copper.mu0; % Wikipedia

iron = ofem_v2.materials.Material;
iron.kappa = 1e7; % Wikipedia: Iron
iron.mu = 1/(5000 * iron.mu0); % Wikipedia: Iron

air = ofem_v2.materials.Material;
air.kappa = 1; % Small conductivity to stabilize the computation realistic would be 1e-12 to 1e-15
air.mu = 1/air.mu0;

mesh.setMaterial('Air', air);
mesh.setMaterial('Coil', copper);
mesh.setMaterial('Core', iron);

%%
upperVoltage = 1;
lowerVoltage = 0;

voltage = ofem_v2.Physical_Problem(fe,mesh,1,0,0);
voltage.setParaS('kappa');
voltage.attachDOFHandler(dofs);

upperBD = ofem_v2.boundary.Dirichlet(upperVoltage,'UpperPort',mesh);
lowerBD = ofem_v2.boundary.Dirichlet(lowerVoltage,'LowerPort',mesh);

voltage.setBoundaryCondition(upperBD);
voltage.setBoundaryCondition(lowerBD);

dofs.generateDOFs;

voltage.assemble;
t = toc;
disp("Done in:")
disp(t)
%voltage.solve;

%%
tic
disp("Solving System")
L = ichol(voltage.S(dofs.freeDOFs,dofs.freeDOFs));

u = full(voltage.u);
u(dofs.freeDOFs) = pcg(voltage.S(dofs.freeDOFs,dofs.freeDOFs),voltage.b(dofs.freeDOFs),1e-9,1000,L,L');
t = toc;
disp("Done in:")
disp(t)
%%
E = -voltage.gradCell(full(u));

kappa = ofem_v2.tools.matrixarray(zeros(1,1,mesh.Nint));
kappa(1,1,mesh.parts{3,1}) = copper.kappa;
kappa(1,1,mesh.parts{3,2}) = iron.kappa;

J = kappa*E;

% mesh.export_UCD([pwd,'/export'],'voltage',{'U',u,''},{'E',squeeze(E)','','Cell'}, {'J',squeeze(J)','','Cell'});

%%
tic;
disp("Setting up and assembling HCurl Problem.")
feNe = ofem_v2.elements.loadFE('HCurl_3D_Order_0');

dofsNe = ofem_v2.DOFHandler(mesh);
dofsNe.attach(feNe);

v_potential = ofem_v2.Physical_Problem(feNe,mesh,1,1,0);
v_potential.attachDOFHandler(dofsNe);

v_potential.setParaS('mu');
v_potential.setParaM('kappa');

upperBD = ofem_v2.boundary.DirichletEdge([0,0,0]','UpperPort',mesh);
lowerBD = ofem_v2.boundary.DirichletEdge([0,0,0]','LowerPort',mesh);

v_potential.setBoundaryCondition(upperBD);
v_potential.setBoundaryCondition(lowerBD);

force = ofem_v2.Volume_force(J,'Coil',mesh);

mesh.setForce('Coil',force);

dofsNe.generateDOFs;
v_potential.assemble();
t = toc;
disp("Done in:")
disp(t)

%%
tic
disp("Solving System.")
S = v_potential.S;
M = v_potential.M;

b = v_potential.b;

omega = 2*pi*1e6; % 1MHz

A = S + 1i*omega*M;
[P,R,C] = equilibrate(A(dofsNe.freeDOFs,dofsNe.freeDOFs));
B = R*P*A(dofsNe.freeDOFs,dofsNe.freeDOFs)*C;
d = R*P*b(dofsNe.freeDOFs);

%L = ichol(M(dofsNe.freeDOFs,dofsNe.freeDOFs));
%L = ichol(B);
%[L,U] = ilu(A(dofsNe.freeDOFs,dofsNe.freeDOFs));
[L,U] = ilu(abs(B));

%%
A_vec = full(v_potential.u);
%gpuB = gpuArray(B);
%gpud = gpuArray(d);
%A_vec(dofsNe.freeDOFs) = bicgstab(A(dofsNe.freeDOFs,dofsNe.freeDOFs),b(dofsNe.freeDOFs),1e-6,20000,L,L');
y = bicgstabl(B,d,1e-6,2000,L,U);
%y = gmres(B,d,300,1e-6,30);
A_vec(dofsNe.freeDOFs) = C*y;
t = toc;
disp("Done in:")
disp(t)

%%

[~,A_rec] = ofem_v2.tools.reconstruct(v_potential,full(A_vec),1);
[~,B] = ofem_v2.tools.reconstructCurl(v_potential,full(A_vec));

%%

H = B;
H(:,mesh.parts{3,1}) = H(:,mesh.parts{3,1})/mesh.parts{2,1}.mu;
H(:,mesh.parts{3,2}) = H(:,mesh.parts{3,2})/mesh.parts{2,2}.mu;
H(:,mesh.parts{3,3}) = H(:,mesh.parts{3,3})/mesh.parts{2,3}.mu;

%%
mesh.export_UCD([pwd,'/export'], 'field', {'B_i',imag(B)','','Cell'},{'A_i',imag(A_rec)','','Cell'},...
										  {'B_r',real(B)','','Cell'},{'A_r',real(A_rec)','','Cell'},...
										  {'B_a',abs(B)','','Cell'},{'A_a',abs(A_rec)','','Cell'},...
										  {'H_a',abs(H)','','Cell'},{'H_i',imag(H)','','Cell'},...
										  {'H_r',real(H)','','Cell'});









