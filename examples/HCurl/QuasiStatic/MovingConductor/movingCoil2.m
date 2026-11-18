%% Copper ring, first-order approximation: source kappa*v x B0 in the ring
% (no velocity term, the background field is not part of the solution)
close all; clear;

res = movingConductor('coilInAir','Conductor','Coil','Kappa',5.8e7, ...
    'Velocity',[10;0;0],'B0',[0;0;1],'Source','firstOrder', ...
    'Solver','gmres','Tol',1e-4);
plotConvergence(res);
exportResults(res,'movingCoil2');
