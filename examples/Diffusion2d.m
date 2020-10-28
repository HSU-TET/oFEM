%% This skript solves a 2D diffusion problem with the ofem.v.2
%The capacitor has two orthogonal plates with a voltage of 5V. Inside the capacitor are 2 different media. 
% left: % right

close all; 
clear all; 

%% Creating the geometry
file = 'Diffusion2d';

mesh= Geometry(); 
mesh.load_from_msh(file);

%% Choosing function space (Element type and order)
element = P1Element(mesh); 

%% Creating Material Classes
paper = Material(); 
paper.my= 200;
paper.c = 12;

textile = Material();
textile.my = 1; 
textile.c = 1;

%% Modeling the Physical Problem
% in this case a parabolic equation with a stiffness and mass matrix with 2
% Dirichlet Boundary Conditions and 2 Neumann Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
dif = Physical_Problem(element, mesh, 1, 1,0);

% assign material to the parts of the geometry
paper.setMaterial(dif, 'Paper');
textile.setMaterial(dif, 'Textile');

% specify which parameter is neccesarry for which Matrix
dif.chooseParameter('Paper',paper.my, 'stiffness');
dif.chooseParameter('Paper', paper.c, 'mass');
dif.chooseParameter('Textile', textile.my, 'stiffness');
dif.chooseParameter('Textile', textile.c, 'mass');

% create Boundary conditions
leftPlate = DirichletBC(0, 'LowCon', mesh); % 0% concentration
rightPlate = DirichletBC(100, 'HighCon', mesh); % 100% concentration

paper2air = NeumannBC(20, 'PaperBoundary', mesh);
textile2air = NeumannBC(9, 'TextileBoundary', mesh);

% assign to the Problem
dif.setBoundaryCondition(leftPlate);
dif.setBoundaryCondition(rightPlate);
dif.setBoundaryCondition(paper2air);
dif.setBoundaryCondition(textile2air);

% set initial conditions 
init = Initial_Condition(0);


%% Assembly of the matrices and solving the equation
dif.assemble(); 


%% Solving and Exporting the Data
dif.parabolic_solve(init, 100, 1); 

