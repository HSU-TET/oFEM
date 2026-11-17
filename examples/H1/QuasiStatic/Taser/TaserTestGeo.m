%% Taser Test Geo
% by Leonard Spiering 07.12.2019 
% 20microsec puls -> frequency of 50kHz 
% PeakVoltage of 500kV
% Material Parameter: dry Skin (thickness?!)
%                       Nerves completely ignored. 
%https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/
close all;
clear all; 

%% Defining Variables
PotentialLeft = 0; %[V]
PotentialRight = 5000; %[V]

%% Importing the Geometry
fileLeg = [pwd,'\geometry\leg3D']; 

mesh = Geometry(); 
mesh.load_from_msh(fileLeg);

%% Choose function Space
element = P1Element(mesh); 

%% Setting up the Model Problem 
% electro quasi steady:
Voltkappa = Physical_Problem(element, mesh, 1,0,0);
Voltepsilon = Physical_Problem(element, mesh, 1,0,0);
%% Boundary Condition
leftElectrode = DirichletBC(PotentialLeft, 'UpperCircle', mesh); 
rightElectrode = DirichletBC(PotentialRight, 'LowerCircle', mesh); 

Voltkappa.setBoundaryCondition(leftElectrode);
Voltkappa.setBoundaryCondition(rightElectrode);
Voltepsilon.setBoundaryCondition(leftElectrode);
Voltepsilon.setBoundaryCondition(rightElectrode);
%% Setting up the Material

e0 = 8.8541878176E-12;
blood = Material();
blood.epsilon = 5.20E+3*e0;
blood.kappa = 7.01E-1; 

muscle = Material(); 
muscle.epsilon =1.01E+4*e0; 
muscle.kappa = 3.52E-1; 

bone = Material(); % Parameter of Cortical bone
bone.epsilon =2.64E+2*e0;
bone.kappa =2.06E-2; 

skin = Material(); 
skin.epsilon = 1.13E+3*e0; 
skin.kappa =2.73E-4; 

fat = Material();
fat.epsilon =1.63E+2*e0; 
fat.kappa =4.33E-2; 

dt = 0.001; 
Voltkappa.chooseParameter('Blood', blood.kappa, 'stiffness');
Voltkappa.chooseParameter('Skin', skin.kappa, 'stiffness');
Voltkappa.chooseParameter('Muscle', muscle.kappa, 'stiffness');
Voltkappa.chooseParameter('Bone', bone.kappa, 'stiffness');
Voltkappa.chooseParameter('Fat', fat.kappa, 'stiffness');

Voltepsilon.chooseParameter('Blood', blood.epsilon/dt, 'stiffness');
Voltepsilon.chooseParameter('Skin', skin.epsilon/dt, 'stiffness');
Voltepsilon.chooseParameter('Muscle', muscle.epsilon/dt, 'stiffness');
Voltepsilon.chooseParameter('Bone', bone.epsilon/dt, 'stiffness');
Voltepsilon.chooseParameter('Fat', fat.epsilon/dt, 'stiffness');

%% Assembly of the matrices
Voltkappa.assemble(); 
Voltepsilon.assemble();

Voltkappa.solve(); 
Voltepsilon.solve(); 

%% Test- Solving as Elliptic Problem
u_init= Initial_Condition(0);
u_i = zeros(mesh.Nco,1);
u_i(Voltkappa.dof) = u_init.value;  %Maximum variable size allowed by the program is exceeded.
n = 100;

A = Voltkappa.S + Voltepsilon.S/dt;
L = ichol(A(Voltkappa.dof,Voltkappa.dof));

for i =1:n+1
    u_imin1 = u_i;
    b_diri = ones(mesh.Nco,1)*pulse(i*dt); % sind das wirklich die dirichlet werte?!
	b_diri(Voltkappa.dof) = 0;
	b_diri = A*b_diri;
    b = -b_diri+Voltepsilon.S/dt * u_imin1;  % ./ !?
    u_i(Voltepsilon.dof) = cgs(A(Voltepsilon.dof,Voltepsilon.dof),b(Voltepsilon.dof),1e-7,1000,L,L');
    
	mesh.export_UCD([pwd, '/exportLeg'],['exportLeg3D',num2str(i-1)], {'u', u_i, ''});
    
end
