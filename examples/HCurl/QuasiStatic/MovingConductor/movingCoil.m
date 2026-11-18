%% Copper ring moving through a uniform field (velocity term)
close all; clear;

res = movingConductor('coilInAir','Conductor','Coil','Kappa',5.8e7, ...
    'Velocity',[1;0;0],'B0',[0;0;1],'Solver','gmres','Tol',1e-4);
plotConvergence(res);
exportResults(res,'movingCoil');
