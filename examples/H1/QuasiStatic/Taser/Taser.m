%% Taser 
% by Leonard Spiering 07.12.2019 
% 20microsec puls -> frequency of 50kHz 
% PeakVoltage of 500kV
% Material Parameter: dry Skin (thickness?!)
%                       Nerves completely ignored. 
%https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/
close all;
clear; 

%% Defining Variables
PotentialLeft = 0; %[V]
PotentialRight = 1; %[V]

%% Importing the Geometry
fileLeg = [pwd,'\geometry\leg3D']; 

mesh = ofem_v2.Geometry(); 
mesh.load_from_msh(fileLeg);
mesh.create_edges;
mesh.create_faces;
mesh.connectFa2Ed;

%% Choose function Space
element = ofem_v2.elements.loadFE('H1_3D_Order_1');
dh = ofem_v2.DOFHandler(mesh);
dh.attach(element);

%% Setting up the Model Problem 
% electro quasi steady:
leg = ofem_v2.Physical_Problem(element, mesh, 1,0,0);
leg.attachDOFHandler(dh);

%% Boundary Condition
leftElectrode = ofem_v2.boundary.Dirichlet(PotentialLeft, 'leftElectrode', mesh); 
rightElectrode = ofem_v2.boundary.Dirichlet(PotentialRight, 'rightElectrode', mesh); 

leg.setBoundaryCondition(leftElectrode);
leg.setBoundaryCondition(rightElectrode);

%% Setting up the Material

dt = 1e-6;
e0 = 8.8541878176E-12;
blood = ofem_v2.materials.Material();
blood.epsilon = 5.20E+3*e0;
blood.kappa = 7.01E-1; 
blood.stiff = blood.kappa + blood.epsilon/dt;

muscle = ofem_v2.materials.Material(); 
muscle.epsilon =1.01E+4*e0; 
muscle.kappa = 3.52E-1;
muscle.stiff = muscle.kappa + muscle.epsilon/dt;

bone = ofem_v2.materials.Material(); % Parameter of Cortical bone
bone.epsilon =2.64E+2*e0;
bone.kappa =2.06E-2;
bone.stiff = bone.kappa + bone.epsilon/dt;

skin = ofem_v2.materials.Material(); 
skin.epsilon = 1.13E+3*e0; 
skin.kappa =2.73E-4;
skin.stiff = skin.kappa + skin.epsilon/dt;

fat = ofem_v2.materials.Material();
fat.epsilon =1.63E+2*e0; 
fat.kappa =4.33E-2;
fat.stiff = fat.kappa + fat.epsilon/dt;
 
% mesh.setMaterial(blood)  // steckt das mat zum part
% Voltkappa.setParaS('kappa');
mesh.setMaterial('Blood',blood);
mesh.setMaterial('Skin',skin);
mesh.setMaterial('Muscle',muscle);
mesh.setMaterial('Bone',bone);
mesh.setMaterial('Fat',fat);

leg.setParaS('stiff');

dh.generateDOFs;

%% Assembly of the matrices
leg.assemble(); 

dofs = dh.freeDOFs;

%% Test- Solving as Elliptic Problem
u_init= ofem_v2.Initial_Condition(0);
% functions ones is not working
u_i = sparse(zeros(mesh.Nco,1));
u_i(dofs) = u_init.value;  %Maximum variable size allowed by the program is exceeded.
n = 100;

%A = Voltkappa.S + Voltepsilon.S/dt; %Kappa halbe?!
%A = leg.S + Voltepsilon.S/dt;
S = leg.S;

L = ichol(S(dofs,dofs));
b_diri = sparse(zeros(mesh.Nco,1));

% Testing if the puls makes problems: 
for i =1:n+1
    u_imin1 = u_i;
    p=pulse(i*dt);
    %b_diri = ones(mesh.Nco,1)*pulse(i*dt);
	b_diri(leftElectrode.nodes) = pulse(i*dt);
	u_i = b_diri;
	b_diri = S*b_diri;
    b = -b_diri+S * u_imin1;  % ./ !?
    u_i(dofs) = cgs(S(dofs,dofs),b(dofs),1e-7,1000,L,L');
    
	mesh.export_UCD([pwd, '/exportLeg'],['exportLeg3D',num2str(i-1)], {'u', u_i, ''});
    
end
% 
% for i =1:n+1
%     i
%     u_imin1 = u_i;
%     p=pulse(i*dt);
%     b_diri = ones(mesh.Nco,1)*pulse(i*dt); % sind das wirklich die dirichlet werte?!
% 	b_diri(Voltkappa.dof) = 0;
% 	b_diri = A*b_diri;
%     b = -b_diri+Voltepsilon.S/dt * u_imin1;  % ./ !?
%     u_i(Voltepsilon.dof) = cgs(A(Voltepsilon.dof,Voltepsilon.dof),b(Voltepsilon.dof),1e-7,1000,L,L');
%     
% 	mesh.export_UCD([pwd, '/exportLeg'],['exportLeg3D',num2str(i-1)], {'u', u_i, ''});
%     
% end


