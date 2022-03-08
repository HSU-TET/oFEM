
%% This skript solves a 3D plate capacitor with the ofem.v.2
%The capacitor has two orthogonal plates with a voltage of 5V. 
%Inside the capacitor are 2 different media. 
% left: % right

close all; 
clear all; 

%% Creating the geometry
file = './geometry/planarCapacitor3D';

mesh= ofem_v2.Geometry(); 
mesh.load_from_msh(file);

mesh.create_edges();
mesh.create_faces();
mesh.connectFa2Ed();

%% Choosing function space (Element type and order)
%element = ofem_v2.P1Element(mesh); 
fe = ofem_v2.elements.loadFE('H1_3D_Order_1');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);

%% Creating Material Classes
paper = ofem_v2.materials.Material(); 
paper.epsilon = 0.0005;

glas = ofem_v2.materials.Material(); 
glas.epsilon = 0.0001; 

%% Modeling the Physical Problem
% in this case an eliptic equation with only a stiffness matrix with 2
% Dirichlet Boundary Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
capacitor = ofem_v2.Physical_Problem(fe, mesh, 1,0,0);
capacitor.attachDOFHandler(dofs);

% assign material to the parts of the geometry
mesh.setMaterial('LeftCap',paper);
mesh.setMaterial('RightCap',glas);

% specify which parameter is neccesarry for which Matrix
capacitor.setParaS('epsilon');

% create Boundary conditions
leftPlate = ofem_v2.boundary.Dirichlet(0, 'LeftDirichlet', mesh); % 0V Potential on the left plate
rightPlate = ofem_v2.boundary.Dirichlet(5, 'RightDirichlet', mesh); % 5V Potential on the right plate

% assign to the Problem
capacitor.setBoundaryCondition(leftPlate);
capacitor.setBoundaryCondition(rightPlate);

%% Assembly of the matrices and solving the equation
dofs.generateDOFs();
capacitor.assemble(); 
capacitor.solve(); 

%% Exporting the Data
mesh.export_UCD([pwd,'/export'],['3DCapacitor'],{'U',capacitor.u,''});



