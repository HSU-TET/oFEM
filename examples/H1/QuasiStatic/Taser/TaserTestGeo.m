%% Taser Test Geo
% by Leonard Spiering 07.12.2019
% 20microsec puls -> frequency of 50kHz
% PeakVoltage of 500kV
% Material Parameter: dry Skin (thickness?!)
%                       Nerves completely ignored.
%https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/
%
% Electro-quasistatic potential in a 3D leg between two ring electrodes
% (UpperCircle, LowerCircle):
%   div((kappa + epsilon d/dt) grad phi) = 0.
% First the purely resistive (kappa) and purely capacitive (epsilon)
% stationary potentials for PotentialLeft/PotentialRight, then the time
% stepping with the upper electrode driven by pulse(t) and the lower one
% grounded. (The original loop set both electrodes to pulse(t), which
% gives a constant potential.)
close all;
clear;

%% Defining Variables
PotentialLeft = 0;     %[V] UpperCircle
PotentialRight = 5000; %[V] LowerCircle
dt = 0.001;            %[s]
n = 100;               % time steps

%% Importing the Geometry
file = fullfile('geometry','leg3D');
if ~isfile([file,'.msh'])
    system(sprintf('gmsh "%s.geo" -3 -v 0 -o "%s.msh"',file,file));
end
mesh = ofem_v2.Geometry();
mesh.load_from_msh(file);
mesh.create_edges;
mesh.create_faces;
mesh.connectFa2Ed;

%% Choose function Space
element = ofem_v2.elements.loadFE('H1_3D_Order_1');
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

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
leg.attachDOFHandler(dofs);
upper = ofem_v2.boundary.Dirichlet(PotentialLeft, 'UpperCircle', mesh);
lower = ofem_v2.boundary.Dirichlet(PotentialRight, 'LowerCircle', mesh);
leg.setBoundaryCondition(upper);
leg.setBoundaryCondition(lower);
dofs.generateDOFs;

%% Assembly and stationary solutions
leg.setParaS('kappa');
leg.assemble(1);
Sk = leg.S;
leg.solve();
uKappa = leg.u;            % purely resistive

leg.setParaS('epsilon');
leg.assemble(1);
Se = leg.S;
leg.solve();
uEpsilon = leg.u;          % purely capacitive

exportDir = fullfile(pwd,'exportLeg');
exportStatic(mesh,exportDir,uKappa,uEpsilon);

%% Time stepping
iUpper = dofs.n2DOF(upper.nodes);
iLower = dofs.n2DOF(lower.nodes);
uD = zeros(dofs.Nd,1);
uD(iLower) = 0;
dirichlet = @(t) setValues(uD,iUpper,pulse(t));

export = @(i,t,u) exportPotential(mesh,exportDir,'exportLeg3D',i,u);
u = dirichlet(0);
export(0,0,u);
u = eqsTimeStepping(Sk,Se,dofs.freeDOFs,dirichlet,u,dt,n,export);
fprintf('TaserTestGeo: %d steps, potential range at the end [%.3g, %.3g] V\n',n,min(u),max(u));

function u = setValues(u,idx,value)
    u(idx) = value;
end

function exportPotential(mesh,folder,name,i,u) %#ok<INUSD> used inside evalc
    evalc('mesh.export_UCD(folder,[name,num2str(i)],{''u'',u,''''})');
end

function exportStatic(mesh,folder,uKappa,uEpsilon) %#ok<INUSD> used inside evalc
    evalc(['mesh.export_UCD(folder,''staticLeg3D'',{''uKappa'',full(uKappa),''V''},' ...
        '{''uEpsilon'',full(uEpsilon),''V''})']);
end
