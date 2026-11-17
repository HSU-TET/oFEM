%% This skript solves a 2D diffusion problem with the ofem.v.2
% Transient diffusion through two layers (paper and textile):
%   c du/dt - div(D grad u) = 0
% with 0% concentration on the left and 100% on the right end and a
% constant flux through the top and bottom surfaces (Neumann).

close all;
clear;

%% Creating the geometry
file = 'Diffusion2d';
if ~isfile([file,'.msh'])
    system(sprintf('gmsh %s.geo -2 -v 0 -o %s.msh',file,file));
end

mesh = ofem_v2.Geometry();
mesh.load_from_msh(file);
mesh.reorderAC;
mesh.create_edges;

%% Choosing function space (Element type and order)
element = ofem_v2.elements.loadFE('H1_2D_Order_1');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

%% Creating Material Classes
% stiff: diffusion coefficient D, mass: capacity c
paper = ofem_v2.materials.Material();
paper.stiff = 200;
paper.mass = 12;

textile = ofem_v2.materials.Material();
textile.stiff = 1;
textile.mass = 1;

%% Modeling the Physical Problem
% in this case a parabolic equation with a stiffness and mass matrix with 2
% Dirichlet Boundary Conditions and 2 Neumann Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
dif = ofem_v2.Physical_Problem(element, mesh, 1, 1,0);
dif.attachDOFHandler(dofs);

% assign material to the parts of the geometry
mesh.setMaterial('Paper',paper);
mesh.setMaterial('Textile',textile);

% specify which parameter is neccesarry for which Matrix
dif.setParaS('stiff');
dif.setParaM('mass');

% create Boundary conditions
leftPlate = ofem_v2.boundary.Dirichlet(0, 'LowCon', mesh); % 0% concentration
rightPlate = ofem_v2.boundary.Dirichlet(100, 'HighCon', mesh); % 100% concentration

feBd = ofem_v2.elements.loadFE('H1_1D_Order_1');
paper2air = ofem_v2.boundary.Neumann(20, feBd, 'PaperBoundary', mesh);
textile2air = ofem_v2.boundary.Neumann(9, feBd, 'TextileBoundary', mesh);

% assign to the Problem
dif.setBoundaryCondition(leftPlate);
dif.setBoundaryCondition(rightPlate);
dif.setBoundaryCondition(paper2air);
dif.setBoundaryCondition(textile2air);

% set initial conditions
init = ofem_v2.Initial_Condition(0);

dofs.generateDOFs;

%% Assembly of the matrices and solving the equation
dif.assemble(1);

%% Solving and Exporting the Data (100 steps with dt = 1)
dif.parabolic_solve(init, 100, 1, '/export', 'time');
