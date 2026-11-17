%% This skript solves a 2D diffusion problem with the ofem.v.2
%The capacitor has two orthogonal plates with a voltage of 5V. Inside the capacitor are 2 different media. 
% left: % right

close all; 
clear all; 

%% Creating the geometry
file = 'Diffusion2d';

mesh = ofem_v2.Geometry(); 
mesh.load_from_msh(file);

%% Choosing function space (Element type and order)
element = ofem_v2.elements.H1Element(2,1);
element.computeBasis;

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

%% Creating Material Classes
paper = ofem_v2.materials.Material(); 
paper.my= 200;
paper.c = 12;

textile = ofem_v2.materials.Material();
textile.my = 1; 
textile.c = 1;

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
dif.setParaS('my');
dif.setParaM('c');

% create Boundary conditions
leftPlate = ofem_v2.boundary.Dirichlet(0, 'LowCon', mesh); % 0% concentration
rightPlate = ofem_v2.boundary.Dirichlet(100, 'HighCon', mesh); % 100% concentration

paper2air = ofem_v2.boundary.Neumann(20, 'PaperBoundary', mesh);
textile2air = ofem_v2.boundary.Neumann(9, 'TextileBoundary', mesh);

% assign to the Problem
dif.setBoundaryCondition(leftPlate);
dif.setBoundaryCondition(rightPlate);
dif.setBoundaryCondition(paper2air);
dif.setBoundaryCondition(textile2air);

% set initial conditions 
init = Initial_Condition(0);

dofs.generateDOFs;


%% Assembly of the matrices and solving the equation
dif.assemble(); 


%% Solving and Exporting the Data
dif.parabolic_solve(init, 100, 1,[pwd,'/export'],['time']); 

