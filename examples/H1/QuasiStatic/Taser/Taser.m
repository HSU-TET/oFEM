%% Taser
% by Leonard Spiering 07.12.2019
% 20microsec puls -> frequency of 50kHz
% PeakVoltage of 500kV
% Material Parameter: dry Skin (thickness?!)
%                       Nerves completely ignored.
%https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/
%
% Electro-quasistatic potential in a 3D leg:
%   div((kappa + epsilon d/dt) grad phi) = 0,
% left electrode driven by pulse(t), right electrode grounded.
close all;
clear;

%% Defining Variables
PotentialRight = 0; %[V] ground
dt = 1e-6;          %[s]
n = 100;            % time steps

%% Importing the Geometry
fileLeg = fullfile(pwd,'geometry','leg3D');
if ~isfile([fileLeg,'.msh'])
    system(sprintf('gmsh "%s.geo" -3 -v 0 -o "%s.msh"',fileLeg,fileLeg));
end

mesh = ofem_v2.Geometry();
mesh.load_from_msh(fileLeg);
mesh.create_edges;
mesh.create_faces;
mesh.connectFa2Ed;

%% Choose function Space
element = ofem_v2.elements.loadFE('H1_3D_Order_1');
dh = ofem_v2.DOFHandler(mesh);
dh.attach(element);

%% Setting up the Material
e0 = ofem_v2.materials.Material.eps0;
blood = ofem_v2.materials.Material();
blood.epsilon = 5.20E+3*e0;
blood.kappa = 7.01E-1;

muscle = ofem_v2.materials.Material();
muscle.epsilon = 1.01E+4*e0;
muscle.kappa = 3.52E-1;

bone = ofem_v2.materials.Material(); % Parameter of Cortical bone
bone.epsilon = 2.64E+2*e0;
bone.kappa = 2.06E-2;

skin = ofem_v2.materials.Material();
skin.epsilon = 1.13E+3*e0;
skin.kappa = 2.73E-4;

fat = ofem_v2.materials.Material();
fat.epsilon = 1.63E+2*e0;
fat.kappa = 4.33E-2;

mesh.setMaterial('Blood',blood);
mesh.setMaterial('Skin',skin);
mesh.setMaterial('Muscle',muscle);
mesh.setMaterial('Bone',bone);
mesh.setMaterial('Fat',fat);

%% Setting up the Model Problem and Boundary Condition
leg = ofem_v2.Physical_Problem(element, mesh, 1,0,0);
leg.attachDOFHandler(dh);
leftElectrode = ofem_v2.boundary.Dirichlet(0, 'leftElectrode', mesh);
rightElectrode = ofem_v2.boundary.Dirichlet(PotentialRight, 'rightElectrode', mesh);
leg.setBoundaryCondition(leftElectrode);
leg.setBoundaryCondition(rightElectrode);
dh.generateDOFs;

%% Assembly of the matrices: conductivity and permittivity
leg.setParaS('kappa');
leg.assemble(1);
Sk = leg.S;
leg.setParaS('epsilon');
leg.assemble(1);
Se = leg.S;

%% Time stepping
left = dh.n2DOF(leftElectrode.nodes);
right = dh.n2DOF(rightElectrode.nodes);
uD = zeros(dh.Nd,1);
uD(right) = PotentialRight;
dirichlet = @(t) setValues(uD,left,pulse(t));

exportDir = fullfile(pwd,'exportLeg');
export = @(i,t,u) exportPotential(mesh,exportDir,'exportLeg3D',i,u);
u_i = dirichlet(0);
export(0,0,u_i);
u_i = eqsTimeStepping(Sk,Se,dh.freeDOFs,dirichlet,u_i,dt,n,export);
fprintf('Taser: %d steps, potential range at the end [%.3g, %.3g] V\n',n,min(u_i),max(u_i));

function u = setValues(u,idx,value)
    u(idx) = value;
end

function exportPotential(mesh,folder,name,i,u) %#ok<INUSD> used inside evalc
    evalc('mesh.export_UCD(folder,[name,num2str(i)],{''u'',u,''''})');
end
