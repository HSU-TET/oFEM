X = mesh.co(1,:,:);
X = double(X(:));
Y = mesh.co(2,:,:);
Y = double(Y(:));
Z = mesh.co(3,:,:);
Z = double(Z(:));

tetramesh(mesh.el(mesh.parts{3,1},:),[X,Y,Z]);

