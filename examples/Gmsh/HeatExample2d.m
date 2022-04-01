%% This skript solves a 2D heat equation with the ofem.v.2
% Model is a steel rod with a set temperature at each end and neuman at the
% sides

close all; 
clear all; 

%% Creating the geometry
file = '.\geometry\HeatExample2d';
%file = 'HeatExample2D';


mesh= ofem_v2.Geometry(); 
mesh.load_from_msh(file);
mesh.reorderAC;
mesh.create_edges;
mesh.create_faces;
mesh.connectFa2Ed;

%% Choosing function space (Element type and order)
%element = ofem_v2.elements.loadFE('H1_2D_Order_3');
element = ofem_v2.elements.H1Element(2,3);
element.computeBasis;
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

%% Creating Material Classes
steel = ofem_v2.materials.Material(); 
steel.lambda= 48;
steel.c = 120;


%% Modeling the Physical Problem
% in this case a parabolic equation with a stiffness and mass matrix with 2
% Dirichlet Boundary Conditions and 2 Neumann Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
heat = ofem_v2.Physical_Problem(element, mesh, 1, 1,0);
heat.attachDOFHandler(dofs);

% assign material to the parts of the geometry
mesh.setMaterial('Steel',steel);

% specify which parameter is neccesarry for which Matrix
heat.setParaS('lambda');
heat.setParaM('c');

% create Boundary conditions
rightPlate = ofem_v2.boundary.Dirichlet(400, 'HighTemp', mesh);

%surface = NeumannBC(11, 'Boundary', mesh);

% assign to the Problem
heat.setBoundaryCondition(rightPlate);
%heat.setBoundaryCondition(surface);

% set initial conditions 
init =  ofem_v2.Initial_Condition(280);

dofs.generateDOFs;

%% Assembly of the matrices and solving the equation
heat.assemble(); 


%% Solving and Exporting the Data
heat.parabolic_solve(init, 100, 1,'/export/heat/','step'); 

