close all
clear

file = 'tempCube';
mesh = ofem_v2.Geometry;
mesh.load_from_msh(file);

fe = ofem_v2.elements.loadFE('H1_3D_Order_1');

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);
dofs.generateDOFs();

square = ofem_v2.materials.Material;
square.stiff = 1;
square.mass = 1;

feBd = ofem_v2.elements.loadFE('H1_2D_Order_1');

% Input in case of heat equation:
% Temperature
% Alpha, a scale factor for the Dirichlet part of the mixed boundary
% Beta, a scale factor for the Neumann part of the boundary
% Element to be used on the boundary, usually the same as the main element
% one dimension lower
% Name of the boundary in the mesh
% The mesh
leftSide = ofem_v2.boundary.Robin(370,1,1,feBd,'Left',mesh);
rightSide = ofem_v2.boundary.Robin(450,1,1,feBd,'Right',mesh);

mesh.setMaterial('Domain',square);

phys = ofem_v2.Physical_Problem(fe,mesh,1,1,0);

phys.setParaS('stiff');
phys.setParaM('mass');

phys.setBoundaryCondition(leftSide);
phys.setBoundaryCondition(rightSide);

phys.attachDOFHandler(dofs);

phys.assemble();

S = phys.S;
M = phys.M;
b = phys.b;
MR = phys.M_robin;

u = phys.u;
u(:) = 200;

mesh.export_UCD([pwd,'/export'],[file,num2str(0)],{'T',u,''});

dt = 1;

LHS = S+M/dt+MR;

for i = 1:100
	RHS = b+M/dt*u;
	
	u = LHS\RHS;
	
	mesh.export_UCD([pwd,'/export'],[file,num2str(i)],{'T',u,''});
end
	
	
	
	
	
	
	
	
