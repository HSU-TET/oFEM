%% Planar capacitor with two dielectrics (paper and glass)
% Electrostatic potential: div(epsilon grad phi) = 0, 10 V on the left and
% 0 V on the right plate.
close all;
clear;

file = 'planarCapacitor';
if ~isfile([file,'.msh'])
    system(sprintf('gmsh %s.geo -2 -v 0 -o %s.msh',file,file));
end

V1 = 10;
V2 = 0;

geom = ofem_v2.Geometry;
geom.load_from_msh(file);
geom.create_edges;

elem = ofem_v2.elements.loadFE('H1_2D_Order_1');
dofs = ofem_v2.DOFHandler(geom);
dofs.attach(elem);

phys = ofem_v2.Physical_Problem(elem,geom,1,0,0);
phys.attachDOFHandler(dofs);

paper = ofem_v2.materials.Material;
paper.epsilon = 2*paper.eps0;

glas = ofem_v2.materials.Material;
glas.epsilon = 8*paper.eps0;

geom.setMaterial('Paper',paper);
geom.setMaterial('Glas',glas);

phys.setParaS('epsilon');

left = ofem_v2.boundary.Dirichlet(V1, 'LeftPlate', geom);
right = ofem_v2.boundary.Dirichlet(V2,'RightPlate',geom);

phys.setBoundaryCondition(left);
phys.setBoundaryCondition(right);

dofs.generateDOFs;
phys.assemble();
phys.solve();

x = double(squeeze(geom.co(1,:,:)));
y = double(squeeze(geom.co(2,:,:)));
trisurf(geom.el,x,y,phys.u)
colorbar();
