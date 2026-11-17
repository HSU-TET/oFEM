%% Taser Simple 2D Forearm
% by Leonard Spiering 07.12.2019
% 20microsec puls -> frequency of 50kHz
% PeakVoltage of 500kV
% Material Parameter: dry Skin (thickness?!)
%                       Nerves completely ignored.
%https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/
%
% Purely capacitive case (kappa = 0) of the electro-quasistatic problem
%   div(epsilon d/dt grad phi) = 0
% in a simplified 2D cross section of homogeneous skin. Two runs: a
% constant 5 V on the left electrode and the voltage pulse(t). The right
% electrode is grounded.
close all;
clear;

%% Defining Variables
PotentialRight = 0; %[V] ground
dt = 1e-8;          %[s]
n = 50;             % time steps

%% Importing the Geometry
file = fullfile('geometry','Simple2D');
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

%% Setting up the Material: skin everywhere
skin = ofem_v2.materials.Material();
skin.epsilon = 1.13E+3*ofem_v2.materials.Material.eps0;
skin.kappa = 2.73E-4;
for i = 1:size(mesh.parts,2)
    mesh.setMaterial(mesh.parts{1,i},skin);
end

%% Setting up the Model Problem and Boundary Condition
arm = ofem_v2.Physical_Problem(element, mesh, 1,0,0);
arm.attachDOFHandler(dofs);
leftElectrode = ofem_v2.boundary.Dirichlet(0, 'leftElectrode', mesh);
rightElectrode = ofem_v2.boundary.Dirichlet(PotentialRight, 'rightElectrode', mesh);
arm.setBoundaryCondition(leftElectrode);
arm.setBoundaryCondition(rightElectrode);
dofs.generateDOFs;

%% Assembly of the matrix (epsilon only)
arm.setParaS('epsilon');
arm.assemble(1);
Se = arm.S;
Sk = sparse(dofs.Nd,dofs.Nd);

left = dofs.n2DOF(leftElectrode.nodes);
right = dofs.n2DOF(rightElectrode.nodes);
uD = zeros(dofs.Nd,1);
uD(right) = PotentialRight;
exportDir = fullfile(pwd,'exportLeg');

%% Constant 5 V on the left electrode
dirichlet = @(t) setValues(uD,left,5);
export = @(i,t,u) exportPotential(mesh,exportDir,'exportSimple2DepsilonConstPulse5',i,u);
u = zeros(dofs.Nd,1);
export(0,0,u);
u = eqsTimeStepping(Sk,Se,dofs.freeDOFs,dirichlet,u,dt,n,export);
fprintf('Simple2D, constant 5 V: potential range [%.3g, %.3g] V\n',min(u),max(u));

%% Voltage pulse on the left electrode
dirichlet = @(t) setValues(uD,left,pulse(t));
export = @(i,t,u) exportPotential(mesh,exportDir,'exportSimple2Depsilon',i,u);
u = zeros(dofs.Nd,1);
export(0,0,u);
u = eqsTimeStepping(Sk,Se,dofs.freeDOFs,dirichlet,u,dt,n,export);
fprintf('Simple2D, pulse: potential range [%.3g, %.3g] V\n',min(u),max(u));

function u = setValues(u,idx,value)
    u(idx) = value;
end

function exportPotential(mesh,folder,name,i,u) %#ok<INUSD> used inside evalc
    evalc('mesh.export_UCD(folder,[name,num2str(i)],{''u'',u,''''})');
end
