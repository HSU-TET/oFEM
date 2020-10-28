X0 = randn(length(dofs.freeDOFs),10);

[S_g,Y_g] = fe.AMS(phys);

% wDofs = TR.freeBoundary;
% wDofs_N0 = wDofs(:);
% wDofs_N0 = unique(wDofs_N0);
% wDofs_E = dofs.boundaryDOFs{1,1};
% wDofs_E = wDofs_E+mesh.Nco;

% wDofs = setdiff(1:mesh.Nco+mesh.Ned,[wDofs_N0;wDofs_E]);

wDofs = dofs.AMSDOFs(phys);

tic
[V,l,Wn,lh] = ofem_v2.solvers.lobpcgModProj(X0,S,M,speye(size(S)),5,300,1e-8,Y_g(:,wDofs),S_g(wDofs,wDofs),dofs.freeDOFs);
toc

v1 = zeros(length(dofs.DOFs),size(V,2));
v1(dofs.freeDOFs,:) = V;


% reconstructTest;