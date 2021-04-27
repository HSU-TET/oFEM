close all;
clear;
inp_file_name='test';


%% set constants
U_1    = 0;
U_2    = 1;

eps1_r = 1; % Left
eps2_r = 2; % Right
eps3_r = 1; % Room


%% load mesh
fprintf('Loading mesh ... \n');
tic
mesh=ofem_v2.Geometry;
mesh.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);


%% define function space discretization
fe = ofem_v2.elements.H1Element(2,1);
fe.computeBasis;

%% create DOFs

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);
dofs.generateDOFs();


leftMat = ofem_v2.materials.Material();
leftMat.epsilon = 5*leftMat.eps0;

rightMat = ofem_v2.materials.Material();
rightMat.epsilon = 5*rightMat.eps0;

roomMat = ofem_v2.materials.Material();
roomMat.epsilon = 5*roomMat.eps0;

leftPlate = ofem_v2.boundary.Dirichlet(1,'LeftPlate',mesh);
rightPlate = ofem_v2.boundary.Dirichlet(0,'RightPlate',mesh);
%testBD = ofem_v2.boundary.Dirichlet(2,'TestBD',mesh);

mesh.setMaterial('LeftSide',leftMat);
mesh.setMaterial('RightSide',rightMat);


phys = ofem_v2.Physical_Problem(fe,mesh,1,0,0);

phys.setParaS('epsilon');

phys.setBoundaryCondition(leftPlate);
phys.setBoundaryCondition(rightPlate);
%phys.setBoundaryCondition(testBD);

phys.attachDOFHandler(dofs);

phys.assemble();

u = phys.u;
b = phys.b;
S = phys.S;

u(dofs.freeDOFs) = S(dofs.freeDOFs,dofs.freeDOFs)\b(dofs.freeDOFs);





% %% define equation type in oFEM 
% op=ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));
% 
% 
% %% oFEM assmble method
% % Use the assemble method of the laplace class to assemble all desired
% % matrizes and vectors.
% fprintf('Assembling stiffness matrix (oFEM) ... ');
% 
% % only stiffness matrix
% opt.S            = 1;
% opt.A            = {eps1_r*mesh.eps0,eps2_r*mesh.eps0,eps3_r*mesh.eps0};
% opt.dirichlet{1} = struct('idx',1,'f',U_1);
% opt.dirichlet{2} = struct('idx',2,'f',U_2);
% [asm,info,aux]     = op.assemble(opt);
% 
% fprintf('done t=%f\n',info.time2assemble);
% 
% 
% %% total assembly and scalar material incorporation
% S    = asm.S;
% DOFs = asm.DOFs;
% b    = asm.b;
% 
% 
% %% solution vector u
% % u = zeros(size(mesh.co,3),1);
% u = full(asm.dirichlet);
% 
% %% solve this guy 
% % NOTE: solve only for DOFs
% fprintf('Solving system ... ');
% 
% tic
% u(DOFs) = S(DOFs,DOFs) \ b(DOFs);
% t=toc;
% fprintf('done t=%f\n',t);
% 
% 
% %% export
mesh.export_UCD(fullfile(pwd,'PlanarCapcitor'),'export',{'U',u,'V'});
% 
% 
% %% plot solution
% figure;
% h=op.plot(u);
% set(h,'FaceColor', 'interp', 'EdgeColor', 'black');

