close all;
clear;

file = 'planarCapacitor';

system(sprintf('gmsh planarCapacitor.geo -setnumber h %d -2 -v 0', 10));

V1 = 10;
V2 = 0;

geom = ofem_v2.Geometry;
geom.load_from_msh(file);

elem = ofem_v2.elements.P1Element(geom);

phys = ofem_v2.Physical_Problem(elem,geom,1,0,0);

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

phys.assemble();
%phys.solve();

S = gpuArray(phys.S(p));
b = gpuArray(phys.b);

x = squeeze(geom.co(1,:,:));
x = double(x);
y = squeeze(geom.co(2,:,:));
y = double(y);
trisurf(geom.el,x,y,phys.u)
colorbar();