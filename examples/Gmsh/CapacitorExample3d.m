
%% This skript solves a 3D plate capacitor with the ofem.v.2
%The capacitor has two orthogonal plates with a voltage of 5V. 
%Inside the capacitor are 2 different media. 
% left: % right

close all; 
clear all; 

%% Creating the geometry
file = '/geometry/planarCapacitor3D';

mesh= ofem_v2.Geometry(); 
mesh.load_from_msh(file);

%% Choosing function space (Element type and order)
element = P1Element(mesh); 

%% Creating Material Classes
paper = Material(); 
paper.epsilon = 0.0005;

glas = Material();
glas.epsilon = 0.0001; 

%% Modeling the Physical Problem
% in this case an eliptic equation with only a stiffness matrix with 2
% Dirichlet Boundary Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
capacitor = Physical_Problem(element, mesh, 1,0,0);

% assign material to the parts of the geometry
paper.setMaterial(capacitor, 'LeftCap');
glas.setMaterial(capacitor, 'RightCap');

% specify which parameter is neccesarry for which Matrix
capacitor.chooseParameter('LeftCap',paper.epsilon, 'stiffness');
capacitor.chooseParameter('RightCap', glas.epsilon, 'stiffness');

% create Boundary conditions
leftPlate = DirichletBC(0, 'LeftDirichlet', mesh); % 0V Potential on the left plate
rightPlate = DirichletBC(5, 'RightDirichlet', mesh); % 5V Potential on the right plate

% assign to the Problem
capacitor.setBoundaryCondition(leftPlate);
capacitor.setBoundaryCondition(rightPlate);

%% Assembly of the matrices and solving the equation
capacitor.assemble(); 
capacitor.solve(); 

%% Exporting the Data
mesh.export_UCD([pwd,'/export'],['exportCap3D'],{'U',capacitor.u,''});



