v12 = mesh.DinvT'*mesh.co(:,1,mesh.el(:,2))-mesh.DinvT'*mesh.co(:,1,mesh.el(:,1));
v12 = v12*(1/ofem_v2.tools.matrixarray(vecnorm(v12)));

v13 = mesh.DinvT'*mesh.co(:,1,mesh.el(:,3))-mesh.DinvT'*mesh.co(:,1,mesh.el(:,1));
v13 = v13*(1/ofem_v2.tools.matrixarray(vecnorm(v13)));

v14 = mesh.DinvT'*mesh.co(:,1,mesh.el(:,4))-mesh.DinvT'*mesh.co(:,1,mesh.el(:,1));
v14 = v14*(1/ofem_v2.tools.matrixarray(vecnorm(v14)));

v23 = mesh.DinvT'*mesh.co(:,1,mesh.el(:,3))-mesh.DinvT'*mesh.co(:,1,mesh.el(:,2));
v23 = v23*(1/ofem_v2.tools.matrixarray(vecnorm(v23)));

v24 = mesh.DinvT'*mesh.co(:,1,mesh.el(:,4))-mesh.DinvT'*mesh.co(:,1,mesh.el(:,2));
v24 = v24*(1/ofem_v2.tools.matrixarray(vecnorm(v24)));

v34 = mesh.DinvT'*mesh.co(:,1,mesh.el(:,3))-mesh.DinvT'*mesh.co(:,1,mesh.el(:,4));
v34 = v34*(1/ofem_v2.tools.matrixarray(vecnorm(v34)));

velem = [squeeze(v12);squeeze(v13);squeeze(v14);squeeze(v23);squeeze(v24);squeeze(v34)];

velem = reshape(velem,3,6,[]);