close all;
clear;

%% input filename
inp_file_name='PlanarCapacitor2D';


%% set constants
U_1    = 1; %leftPlate
U_2    = 0; %rightPlate

eps1_r = 1; % Left
eps2_r = 2; % Right


%% load mesh
fprintf('Loading mesh ... \n');
tic
mesh=ofem_v2.Geometry;
mesh.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);


%% define function space discretization
fe = ofem_v2.elements.loadFE('H1_2D_Order_1');
%fe.computeBasis;


%% create DOFs
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);
dofs.generateDOFs();


%% set materials
leftMat = ofem_v2.materials.Material();
mateqleft = @(X) eps1_r*leftMat.eps0;     %*abs(double(X(1,:,:)));
leftMat.epsilon = mateqleft;

rightMat = ofem_v2.materials.Material();
mateqright = @(X) eps2_r*leftMat.eps0;    %*abs(double(X(1,:,:)));
rightMat.epsilon = mateqright;

mesh.setMaterial('LeftSide',leftMat);
mesh.setMaterial('RightSide',rightMat);


%% define boundary conditions (function handle or constant)
% leftPlate = ofem_v2.boundary.Dirichlet(@leftPlate,'LeftPlate',mesh);
leftPlate = ofem_v2.boundary.Dirichlet(U_1,'LeftPlate',mesh);
rightPlate = ofem_v2.boundary.Dirichlet(U_2,'RightPlate',mesh);


%% oFEM assemble method
phys = ofem_v2.Physical_Problem(fe,mesh,1,0,0);

phys.setParaS('epsilon');

%set boundary conditions
phys.setBoundaryCondition(leftPlate);
phys.setBoundaryCondition(rightPlate);

%set DOFs
phys.attachDOFHandler(dofs);

%assemble
phys.assemble();

%save needed vectors and matrices
u = phys.u;
b = phys.b;
S = phys.S;

%solve
u(dofs.freeDOFs) = S(dofs.freeDOFs,dofs.freeDOFs)\b(dofs.freeDOFs);

%% export
mesh.export_UCD(fullfile(pwd,'PlanarCapcitor2D'),'export',{'U',full(u),'V'});
 
 
%% plot solution (only 2D)
figure;
h=ofem_v2.tools.plot(mesh,u,'V');
set(h,'FaceColor', 'interp', 'EdgeColor', 'black');