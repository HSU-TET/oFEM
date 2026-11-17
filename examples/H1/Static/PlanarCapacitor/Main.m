%% Poisson problem with Dirichlet and Neumann data and a volume force
% -div(rho grad u) = f in the planar capacitor geometry with
%   u = 10            on the left plate   (Dirichlet)
%   rho du/dn = 30    on the right plate  (Neumann)
%   f = 9             in both regions     (volume force)
% and rho = 2 in the paper, rho = 100 in the glass region. (Ported from the
% first oFEM version, where the same setup was used on a wire mesh.)
close all;
clear;

%% Setting up the Geometry class
file = 'planarCapacitor';
if ~isfile([file,'.msh'])
    system(sprintf('gmsh %s.geo -2 -v 0 -o %s.msh',file,file));
end
mesh = ofem_v2.Geometry();
mesh.load_from_msh(file);
mesh.create_edges;

%% Choosing the function Space
element = ofem_v2.elements.loadFE('H1_2D_Order_1');
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(element);

%% Creating Material Classes
surrounding = ofem_v2.materials.Material();
surrounding.rho = 2;
wire = ofem_v2.materials.Material();
wire.rho = 100;
mesh.setMaterial('Paper', surrounding);
mesh.setMaterial('Glas', wire);

%% Modeling the Physical Problem
problem = ofem_v2.Physical_Problem(element, mesh, 1,0,0);
problem.attachDOFHandler(dofs);
problem.setParaS('rho');

%% Boundary conditions and volume force
feBd = ofem_v2.elements.loadFE('H1_1D_Order_1');
dirichlet = ofem_v2.boundary.Dirichlet(10, 'LeftPlate', mesh);
neumann = ofem_v2.boundary.Neumann(30, feBd, 'RightPlate', mesh);
problem.setBoundaryCondition(dirichlet);
problem.setBoundaryCondition(neumann);
ofem_v2.Volume_force(9, 'Paper', mesh);
ofem_v2.Volume_force(9, 'Glas', mesh);

%% Assembly and solution
dofs.generateDOFs;
problem.assemble();
problem.solve();

%% Exporting the Data
mesh.export_UCD([pwd,'/export'],'export',{'U',problem.u,''});
