%% Taser Simple 2D Forearm
% by Leonard Spiering 07.12.2019
% 20microsec puls -> frequency of 50kHz
% PeakVoltage of 500kV
% Material Parameter: dry Skin (thickness?!)
%                       Nerves completely ignored.
%https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/
%
% Electro-quasistatic potential in a simplified 2D cross section (two
% layers), left electrode driven by pulse(t), right electrode grounded.
% Exports the potential and the cell-wise field E = -grad(phi) and
% current density J = kappa*E.
close all;
clear;

%% Defining Variables
PotentialRight = 0; %[V] ground
dt = 1e-3;          %[s]
n = 500;            % time steps

%% Importing the Geometry
file = fullfile('geometry','Simple2D');
if ~isfile([file,'.msh'])
    system(sprintf('gmsh "%s.geo" -2 -v 0 -o "%s.msh"',file,file));
end
mesh = ofem_v2.Geometry();
mesh.load_from_msh(file);
mesh.create_edges;
mesh.jacobiandata;

%% Choose function Space
element = ofem_v2.elements.loadFE('H1_2D_Order_1');
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

%% Setting up the Material
e0 = ofem_v2.materials.Material.eps0;
muscle = ofem_v2.materials.Material();
muscle.epsilon = 1.01E+4*e0;
muscle.kappa = 3.52E-1;

fat = ofem_v2.materials.Material();
fat.epsilon = 1.63E+2*e0;
fat.kappa = 4.33E-2;

% as in the original test: fat properties in the outer layer, muscle
% properties in the inner one
mesh.setMaterial('skin',fat);
mesh.setMaterial('fat',muscle);

%% Setting up the Model Problem and Boundary Condition
arm = ofem_v2.Physical_Problem(element, mesh, 1,0,0);
arm.attachDOFHandler(dofs);
leftElectrode = ofem_v2.boundary.Dirichlet(0, 'leftElectrode', mesh);
rightElectrode = ofem_v2.boundary.Dirichlet(PotentialRight, 'rightElectrode', mesh);
arm.setBoundaryCondition(leftElectrode);
arm.setBoundaryCondition(rightElectrode);
dofs.generateDOFs;

%% Assembly of the matrices
arm.setParaS('kappa');
arm.assemble(1);
Sk = arm.S;
arm.setParaS('epsilon');
arm.assemble(1);
Se = arm.S;

% conductivity per triangle for J = kappa*E
kappa = zeros(1,1,mesh.Nint);
for i = 1:size(mesh.parts,2)
    kappa(mesh.parts{3,i}) = mesh.parts{2,i}.kappa;
end

%% Time stepping
left = dofs.n2DOF(leftElectrode.nodes);
right = dofs.n2DOF(rightElectrode.nodes);
uD = zeros(dofs.Nd,1);
uD(right) = PotentialRight;
dirichlet = @(t) setValues(uD,left,pulse(t));

exportDir = fullfile(pwd,'exportTest');
export = @(i,t,u) exportFields(mesh,arm,kappa,exportDir,i,u);
u = dirichlet(0);
export(0,0,u);
u = eqsTimeStepping(Sk,Se,dofs.freeDOFs,dirichlet,u,dt,n,export);
fprintf('simpleTest: %d steps, potential range at the end [%.3g, %.3g] V\n',n,min(u),max(u));

function u = setValues(u,idx,value)
    u(idx) = value;
end

function exportFields(mesh,phys,kappa,folder,i,u) %#ok<INUSL> mesh, folder, i used inside evalc
    E = -phys.gradCell(u);                  % 2 x 1 x Ne
    J = kappa.*E;
    Ne = size(E,3);
    Eexp = [reshape(E,2,[])', zeros(Ne,1)]; %#ok<NASGU> used inside evalc
    Jexp = [reshape(J,2,[])', zeros(Ne,1)]; %#ok<NASGU> used inside evalc
    evalc(['mesh.export_UCD(folder,[''exportCGS'',num2str(i)],{''u'',u,''''},' ...
        '{''E'',Eexp,'''',''Cell''},{''J'',Jexp,'''',''Cell''})']);
end
