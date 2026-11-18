%% Nickel steel ring moving through the earth's field
close all; clear;

Bdir = [0;30;-60];
B0 = 50e-6*Bdir/norm(Bdir);             % earth's field

res = movingConductor('coilInAir','Conductor','Coil','Mur',40,'Kappa',1e6, ...
    'Velocity',[10;0;0],'B0',B0,'Tol',1e-4);
fprintf('assembly %.1fs, solve %.1fs\n',res.times.assembly,res.times.solve);
plotConvergence(res);
exportResults(res,'movingCoilSteel');
