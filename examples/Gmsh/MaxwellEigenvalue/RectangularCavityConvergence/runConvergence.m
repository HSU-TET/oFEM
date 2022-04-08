close all;
clear;

nh = 10;
no = 3;
err = zeros(20,nh,no+1);

file = 'rect';

for i = 1:nh
	system(sprintf('gmsh.exe rect.geo -setnumber h %d -3 -v 0', i));
	
	mesh = ofem_v2.Geometry;
	mesh.load_from_msh(file);
	
	mesh.reorderAC;
	mesh.create_faces;
	mesh.create_edges;
	mesh.connectFa2Ed;
	mesh.jacobiandata;
	air = ofem_v2.materials.Material;
	air.stiff = 1;
	air.mass = 1;
	if i == 4
		no = no-1;
	end
	if i == 5
		no = no-1;
	end

	
	mesh.setMaterial('Domain',air);
	for j = 0:no
		[err(:,i,j+1),n(i,j+1)] = rect(mesh,j);
	end
end

%%
figure;
err = abs(err);
loglog(n(:,1),err(2,:,1))
hold on
loglog(n(:,2),err(2,:,2))
loglog(n(1:4,3),err(2,1:4,3))
loglog(n(1:3,4),err(2,1:3,4))
legend('0','1','2','3')