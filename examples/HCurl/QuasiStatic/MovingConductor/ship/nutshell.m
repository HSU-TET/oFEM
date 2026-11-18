%% Nutshell hull moving through a uniform field
close all; clear;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % movingConductor

v = [0;1e-4;0];
res = movingConductor('nutshell','Conductor','Hull','Mur',40,'Kappa',1e6, ...
    'Velocity',v,'B0',[0;0;1],'MaxIt',400);
fprintf('assembly %.1fs, solve %.1fs\n',res.times.assembly,res.times.solve);
plotConvergence(res,sprintf('Convergence at $|v|=10^{%d}$',round(log10(norm(v)))));
exportResults(res,'nutshell');
