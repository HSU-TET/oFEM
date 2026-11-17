%% Taser 2D Forearm
% by Leonard Spiering 07.12.2019
% 20microsec puls -> frequency of 50kHz
% PeakVoltage of 500kV
% Material Parameter: dry Skin (thickness?!)
%                       Nerves completely ignored.
%https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/
%
% Electro-quasistatic potential in a cross section of a forearm:
%   div((kappa + epsilon d/dt) grad phi) = 0,
% left electrode driven by pulse(t), right electrode grounded.
close all;
clear;

%% Defining Variables
PotentialRight = 0; %[V] ground
dt = 1e-8;          %[s]
n = 100;            % time steps

%% Importing the Geometry
file = fullfile('geometry','forearm2D');
if ~isfile([file,'.msh'])
    system(sprintf('gmsh "%s.geo" -2 -v 0 -o "%s.msh"',file,file));
end
mesh = ofem_v2.Geometry();
mesh.load_from_msh(file);
mesh.create_edges;

%% Choose function Space
element = ofem_v2.elements.loadFE('H1_2D_Order_1');
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

%% Setting up the Material
e0 = ofem_v2.materials.Material.eps0;

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

mesh.setMaterial('skin',skin);
mesh.setMaterial('muscle',muscle);
mesh.setMaterial('bone',bone);
mesh.setMaterial('fat',fat);

%% Setting up the Model Problem
% electro quasi steady: stiffness matrices with kappa and with epsilon
forearm = ofem_v2.Physical_Problem(element, mesh, 1,0,0);
forearm.attachDOFHandler(dofs);

%% Boundary Condition
leftElectrode = ofem_v2.boundary.Dirichlet(0, 'leftElectrode', mesh);
rightElectrode = ofem_v2.boundary.Dirichlet(PotentialRight, 'rightElectrode', mesh);
forearm.setBoundaryCondition(leftElectrode);
forearm.setBoundaryCondition(rightElectrode);
dofs.generateDOFs;

%% Assembly of the matrices
forearm.setParaS('kappa');
forearm.assemble(1);
Sk = forearm.S;
forearm.setParaS('epsilon');
forearm.assemble(1);
Se = forearm.S;

%% Time stepping
left = dofs.n2DOF(leftElectrode.nodes);
right = dofs.n2DOF(rightElectrode.nodes);
uD = zeros(dofs.Nd,1);
uD(right) = PotentialRight;
dirichlet = @(t) setValues(uD,left,pulse(t));

exportDir = fullfile(pwd,'exportLeg');
export = @(i,t,u) exportPotential(mesh,exportDir,'exportForearm2D',i,u);
u = dirichlet(0);
export(0,0,u);
u = eqsTimeStepping(Sk,Se,dofs.freeDOFs,dirichlet,u,dt,n,export);
fprintf('Taser2D: %d steps, max potential at the end %.3e V\n',n,max(abs(u)));

function u = setValues(u,idx,value)
    u(idx) = value;
end

function exportPotential(mesh,folder,name,i,u) %#ok<INUSD> used inside evalc
    evalc('mesh.export_UCD(folder,[name,num2str(i)],{''u'',u,''''})');
end
