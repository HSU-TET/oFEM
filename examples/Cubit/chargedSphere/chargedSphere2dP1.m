close all;
clear all;

inp_file_name='chargedSphere'; 

%% load mesh
fprintf('Loading mesh ... ');
tic
mesh=ofem_v2.Geometry;
mesh.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);

% figure
% s=mesh.show();
% set(s,'FaceAlpha',0.1);


%% define function space discretization
fe = ofem_v2.elements.H1Element(3,1);
fe.computeBasis;


%% create DOFs
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);
dofs.generateDOFs();

%% define boundary conditions (function handle or constant)
ground = ofem_v2.boundary.Dirichlet(0,'Erde',mesh);

%% set materials
outside= ofem_v2.materials.Material();
outside.epsilon = 1;

sphere = ofem_v2.materials.Material(); 
sphere.epsilon = 1;

mesh.setMaterial('Sphere_1',sphere);
mesh.setMaterial('Aussenraum',outside);

%% Force
ofem_v2.Volume_force(@chargeDensity,1,mesh)


%% oFEM assemble method
phys = ofem_v2.Physical_Problem(fe,mesh,1,0,0);

phys.setParaS('epsilon');

%set boundary conditions
phys.setBoundaryCondition(ground);

%set DOFs
phys.attachDOFHandler(dofs);

%assemble
phys.assemble();

%save needed vectors and matrices
u = phys.u;
b = phys.b;
S = phys.S;

%solve
u(dofs.freeDOFs) = S(dofs.freeDOFs,dofs.freeDOFs)\b(dofs.freeDOFs);

%% compute gradient
E= -phys.gradu(full(u));

 %% export
mesh.export_UCD(fullfile(pwd,'chargedSphere'), 'sphere', {'U', full(u), 'V'}, ...
                {'E', E,'V/m' });
 
%% plot solution

inner_radius = norm(mesh.co(:,1,mesh.bd{2,2}(1,1)));
outer_radius = norm(mesh.co(:,1,mesh.bd{2,1}(1,1)));

% Plot the E_Abs
figure;
plot(squeeze(sqrt(dot(mesh.co,mesh.co,1))),sqrt(dot(E,E,2)),'m+');
legend('E_{Abs}')
xlabel('Distance')
title(['Dimensions: inner radius = ',num2str(inner_radius),', outer radius = ',num2str(outer_radius)]);