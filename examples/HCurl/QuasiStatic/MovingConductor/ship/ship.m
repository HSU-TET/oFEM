%% Steel hull moving through a uniform field
close all; clear;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % movingConductor

res = movingConductor('ship','Conductor','Hull','Mur',40,'Kappa',1e6, ...
    'Velocity',[0;5;0],'B0',[0;0;1]);
plotConvergence(res);
exportResults(res,'ship');
