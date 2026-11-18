%% Uniform background field in air (test of the IABC and the field source)
% The field inside the IABC shells has to be B0 everywhere.
close all; clear;

B0 = [0;0;1];
res = movingConductor('airOnly','Conductor','','B0',B0,'Solver','gmres','Tol',1e-6);

err = vecnorm(res.Bind)/norm(B0);
fprintf('IABC shells mu_r = %s (R = %g, d = %g)\n',mat2str(res.mur',6),res.R,res.d);
fprintf('|B - B0|/|B0| in the air: max %.2e, mean %.2e\n', ...
    max(err(res.mesh.parts{3,strcmp(res.mesh.parts(1,:),'Air')})), ...
    mean(err(res.mesh.parts{3,strcmp(res.mesh.parts(1,:),'Air')})));
exportResults(res,'airOnly');
