% Generate a simple Test Mesh, to check geometry functionality
close all
% Instantiate geometry
mesh = ofem_v2.Geometry;

% Fill basic data
mesh.dim = 3;
mesh.type = 'tet';

%% Use a [0,1]^3 Cube, here with 6 tets
%x = 0:0.5:1;
%y = 0:0.5:1;
%z = 0:0.5:1;

%% Switch to cavity Jin page 272
% x = linspace(0,0.01,11);
% y = linspace(0,0.005,7);
% z = linspace(0,0.0075,11);

disc = 10;

x = linspace(0,pi,disc);
y = linspace(0,pi,disc);
z = linspace(0,pi,disc);

[X,Y,Z] = meshgrid(x,y,z);

X = X(:);
Y = Y(:);
Z = Z(:);

TR = delaunayTriangulation(X,Y,Z);

co = ofem_v2.tools.matrixarray(reshape([X,Y,Z]',3,1,[]));

mesh.el = TR.ConnectivityList;
mesh.co = ofem_v2.tools.matrixarray(reshape(TR.Points',3,1,[]));

mesh.Nco = size(mesh.co,3);

mesh.Nint = size(mesh.el,1);

%clear co el x X y Y z Z

%% Next fill the rest of the structure

mesh.reorderAC;

mesh.create_edges;
mesh.create_faces;

mesh.connectFa2Ed;

mesh.parts{1,1} = 'Cube';
mesh.parts{3,1} = 1:size(mesh.el,1);
mesh.bd{1,1} = 'Left';
mesh.bd{2,1} = TR.freeBoundary;

mesh.jacobiandata;

%% Generate the element
% fe = ofem_v2.elements.NedelecZagl(3,4);
% [N,curlN] = fe.computeBasis;
fe = ofem_v2.elements.loadFE('NE3');

%% Generate DOFs
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);
dofs.generateDOFs();

%% Materials
mat = ofem_v2.materials.Material();
mat.kappa = 1;

%% Boundaries
PEC = ofem_v2.boundary.DirichletEdge([1;0;0],'Left',mesh);

%% Start assembling
phys = ofem_v2.Physical_Problem(fe,mesh,1,1,0);
mesh.setMaterial('Cube',mat);
phys.setParaS('kappa');
phys.setParaM('kappa');
phys.setBoundaryCondition(PEC);
phys.attachDOFHandler(dofs);

%%
phys.assemble();

u = phys.u;
% reconstructTest;

S = phys.S;%(dofs.freeDOFs,dofs.freeDOFs);
M = phys.M;%(dofs.freeDOFs,dofs.freeDOFs);
% 
% opts.p = 20;
% opt.tol = 1e-14;
% 
% tic
% [v,temp] = eigs(S(dofs.freeDOFs,dofs.freeDOFs),M(dofs.freeDOFs,dofs.freeDOFs),15,2,opts);
% toc
% ksq = diag(temp);
% % ksq(:,i) = diag(temp);
% % i=i+1;
% % %k2 = eigs(S,M,20,-1e7);
% v1 = zeros(length(dofs.DOFs),size(v,2));
% v1(dofs.freeDOFs,:) = v;
% 
% 
% %reconstructTest;


femLOBPCG






