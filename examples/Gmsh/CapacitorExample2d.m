
%% This skript solves a 2D plate capacitor with the ofem.v.2
%The capacitor has two orthogonal plates with a voltage of 5V. 
%Inside the capacitor are 2 different media. 
% left: % right

close all; 
clear all; 

%% Creating the geometry
file = './geometry/CapacitorExample2D';

mesh= ofem_v2.Geometry(); 
mesh.load_from_msh(file);
mesh.create_edges;

%% Choosing function space (Element type and order)
%element = ofem_v2.elements.H1Element(2,1); 
element = ofem_v2.elements.loadFE('H1_2D_Order_1');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

%% Creating Material Classes
paper = ofem_v2.materials.Material(); 
paper.epsilon = 0.0001;

glas = ofem_v2.materials.Material();
glas.epsilon = 0.0005; 

%% Modeling the Physical Problem
% in this case an eliptic equation with only a stiffness matrix with 2
% Dirichlet Boundary Conditions

%objName = Constructor(Element, Geometry, has_Stiffness, has_Mass, has_Damping)
capacitor = ofem_v2.Physical_Problem(element, mesh, 1,0,0);
capacitor.attachDOFHandler(dofs)

% assign material to the parts of the geometry
%paper.setMaterial(capacitor, 'Paper');
%glas.setMaterial(capacitor, 'Glas');

mesh.setMaterial('Paper', paper);
mesh.setMaterial('Glas', glas)

% specify which parameter is neccesarry for which Matrix
%capacitor.chooseParameter('Paper',paper.epsilon, 'stiffness');
%capacitor.chooseParameter('Glas', glas.epsilon, 'stiffness');
capacitor.setParaS('epsilon');


% create Boundary conditions
leftPlate = ofem_v2.boundary.Dirichlet(0, 'left', mesh); % 0V Potential on the left plate
rightPlate = ofem_v2.boundary.Dirichlet(5, 'right', mesh); % 5V Potential on the right plate

% assign to the Problem
capacitor.setBoundaryCondition(leftPlate);
capacitor.setBoundaryCondition(rightPlate);

dofs.generateDOFs;

%% Assembly of the matrices and solving the equation
capacitor.assemble(); 
capacitor.solve(); 

%% Exporting the Data
mesh.export_UCD([pwd,'/export'],['2dCapacitor'],{'U',capacitor.u,''});



