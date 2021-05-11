
close all; 
clear all; 

%% Setting up the Geometry class
file = 'planarCapacitor';  % 2D mesh from exercise 1
% load file
mesh= Geometry(); 
mesh.load_from_msh(file);
% Setting up Partname and BD names
mesh.parts{1,1}= 'Wire';
mesh.parts{1,2}= 'Surrounding';
mesh.bd{1,1}= 'BD:1';
mesh.bd{1,2}= 'BD:2';

%% Choosing the function Space

%creating object of P1Element class
element = P1Element(mesh);


%% Creating Material Classes

%create matetial class
air = Material(); 
cupper =Material();
%assign values to the properties
air.rho = 2; 
air.kappa = 4; 
air.my = 7;
cupper.rho = 100; 
cupper.kappa = 400; 
cupper.my= 700;

%% Modeling the Physical Problem

%creating class with functionspace, mesh and specify matrices
Cap2d = Physical_Problem(element, mesh, 1,0,0); 
% assign material to the parts
air.setMaterial(Cap2d, 'Surrounding');
cupper.setMaterial(Cap2d, 'Wire');

%choose which Materialparameter is needed for which Matrix
Cap2d.chooseParameter('Surrounding',air.rho, 'stiffness');
Cap2d.chooseParameter('Surrounding',air.kappa, 'mass');
Cap2d.chooseParameter('Surrounding',air.my, 'damping');

Cap2d.chooseParameter('Wire',cupper.rho, 'stiffness');
Cap2d.chooseParameter('Wire',cupper.kappa, 'mass');
Cap2d.chooseParameter('Wire',cupper.my, 'damping');
%% Creating BC

dirichlet_1 = DirichletBC(10, 'BD:1', mesh);
neumann_1 = NeumannBC(30, 'BD:2', mesh);
%dirichlet_2 = DirichletBC(15, 'BD:2', mesh);

%% Setting Volume Force and BC
Cap2d.volume_force= Volume_force(9, [1,1,2,3,4,5]); 

Cap2d.setBoundaryCondition(dirichlet_1);
Cap2d.setBoundaryCondition(neumann_1);

%% Assembly of the matrices

Cap2d.assemble(); 

%% Solving the System
Cap2d.solve(); 

%% Exporting the Data
mesh.export_UCD([pwd,'/export'],['export'],{'U',Cap2d.u,''});






