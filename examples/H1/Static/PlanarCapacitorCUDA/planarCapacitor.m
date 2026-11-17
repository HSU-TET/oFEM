%% Planar capacitor solved on the GPU
% Same problem as ../PlanarCapacitor: div(epsilon grad phi) = 0 with two
% dielectrics, 10 V on the left and 0 V on the right plate. The reduced
% system is solved with pcg on the GPU (Parallel Computing Toolbox) if
% one is available, otherwise on the CPU.
close all;
clear;

file = 'planarCapacitor';
h = 10;     % mesh size height/h, increase it for a problem worth the GPU
system(sprintf('gmsh %s.geo -setnumber h %d -2 -v 0 -o %s.msh',file,h,file));

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

%% Solve the reduced system with pcg on the GPU
free = dofs.freeDOFs;
S = phys.S(free,free);
b = full(phys.b(free));
s = 1/max(abs(diag(S)));          % scale, epsilon is ~1e-11
if canUseGPU
    S = gpuArray(S*s);
    b = gpuArray(b*s);
    device = 'GPU';
else
    S = S*s;
    b = b*s;
    device = 'CPU (no GPU available)';
end
[x,flag,relres,iter] = pcg(S,b,1e-10,5000);
phys.u(free) = gather(x);
fprintf('pcg on the %s: flag %d, relative residual %.1e, %d iterations\n',device,flag,relres,iter);

x = double(squeeze(geom.co(1,:,:)));
y = double(squeeze(geom.co(2,:,:)));
trisurf(geom.el,x,y,phys.u)
colorbar();
