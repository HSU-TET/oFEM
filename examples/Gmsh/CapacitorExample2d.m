
%% This skript solves a 2D plate capacitor with the ofem.v.2
%The capacitor has two orthogonal plates with a voltage of 5V. 
%Inside the capacitor are 2 different media. 
% left: % right

close all; 
clear all; 

%% Creating the geometry
file = 'CapacitorExample2d';

mesh= Geometry(); 
mesh.load_from_msh(file);

%% Choosing function space (Element type and order)
element = P1Element(mesh); 

%% Creating Material Classes
paper = Material(); 
paper.epsilon = 0.0001;

glas = Material();
glas.epsilon = 0.0005; 

%% Modeling the Physical Problem
% in this case an eliptic equation with only a stiffness matrix with 2
% Dirichlet Boundary Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
capacitor = Physical_Problem(element, mesh, 1,0,0);

% assign material to the parts of the geometry
paper.setMaterial(capacitor, 'Paper');
glas.setMaterial(capacitor, 'Glas');

% specify which parameter is neccesarry for which Matrix
capacitor.chooseParameter('Paper',paper.epsilon, 'stiffness');
capacitor.chooseParameter('Glas', glas.epsilon, 'stiffness');

% create Boundary conditions
leftPlate = DirichletBC(0, 'left', mesh); % 0V Potential on the left plate
rightPlate = DirichletBC(5, 'right', mesh); % 5V Potential on the right plate

% assign to the Problem
capacitor.setBoundaryCondition(leftPlate);
capacitor.setBoundaryCondition(rightPlate);

%% Assembly of the matrices and solving the equation
capacitor.assemble(); 
capacitor.solve(); 

%% Exporting the Data
mesh.export_UCD([pwd,'/export'],['export'],{'U',capacitor.u,''});



