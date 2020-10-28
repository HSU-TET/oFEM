%% This skript solves a 2D heat equation with the ofem.v.2
% Model is a steel rod with a set temperature at each end and neuman at the
% sides

close all; 
clear all; 

%% Creating the geometry
file = 'D:\GitHub\oFEM_rework\geometry\planarCapacitor3D';
%file = 'HeatExample2D';


mesh= Geometry(); 
mesh.load_from_msh(file);

%% Choosing function space (Element type and order)
element = P1Element(mesh); 

%% Creating Material Classes
steel = Material(); 
steel.lambda= 48;
steel.c = 120;


%% Modeling the Physical Problem
% in this case a parabolic equation with a stiffness and mass matrix with 2
% Dirichlet Boundary Conditions and 2 Neumann Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
heat = Physical_Problem(element, mesh, 1, 1,0);

% assign material to the parts of the geometry
steel.setMaterial(heat, 'Steel');


% specify which parameter is neccesarry for which Matrix
heat.chooseParameter('Steel',steel.lambda, 'stiffness');
heat.chooseParameter('Steel', steel.c, 'mass');


% create Boundary conditions
rightPlate = DirichletBC(400, 'HighTemp', mesh);

%surface = NeumannBC(11, 'Boundary', mesh);

% assign to the Problem
heat.setBoundaryCondition(rightPlate);
%heat.setBoundaryCondition(surface);

% set initial conditions 
init = Initial_Condition(280);


%% Assembly of the matrices and solving the equation
heat.assemble(); 


%% Solving and Exporting the Data
heat.parabolic_solve(init, 100, 1); 

