%[N,curlN] = fe.computeBasis();

[X,Y,Z] = meshgrid(0:0.1:1);

% f = matlabFunction(N(1,1));
% g = matlabFunction(N(2,1));
% h = matlabFunction(N(3,1));
f = matlabFunction(N(:,24));
g = matlabFunction(N(:,27));

X(X+Y+Z>1) = NaN;
Y(X+Y+Z>1) = NaN;
Z(X+Y+Z>1) = NaN;

DinvT = mesh.DinvT;

el1 = 4;
el2 = 5;

V = -f(X(:)',Y(:)',Z(:)');
Vnew = -g(X(:)',Y(:)',Z(:)');

V1 = DinvT(:,:,el1)*V;
V2 = DinvT(:,:,el2)*Vnew;
Xnew = double(mesh.co(:,:,mesh.el(el1,1)))+double(mesh.Dk(:,:,el1))*[X(:)';Y(:)';Z(:)'];
X2 = double(mesh.co(:,:,mesh.el(el2,1)))+double(mesh.Dk(:,:,el2))*[X(:)';Y(:)';Z(:)'];

V1 = reshape(V1,3*size(X,1),size(X,2),size(X,3));
V2 = reshape(V2,3*size(X,1),size(X,2),size(X,3));
Xnew = reshape(Xnew,3*size(X,1),size(X,2),size(X,3));
X2 = reshape(X2,3*size(X,1),size(X,2),size(X,3));

% X = X(:,:,1);
% Y = Y(:,:,1);
% Z = Z(:,:,1);
figure
quiver3(Xnew(1:3:end,:,:),Xnew(2:3:end,:,:),Xnew(3:3:end,:,:),V1(1:3:end,:,:),V1(2:3:end,:,:),V1(3:3:end,:,:),1,'LineWidth',2);
hold on
quiver3(X2(1:3:end,:,:),X2(2:3:end,:,:),X2(3:3:end,:,:),V2(1:3:end,:,:),V2(2:3:end,:,:),V2(3:3:end,:,:),1,'LineWidth',2);

V = sqrt(f(X,Y,Z).^2+g(X,Y,Z).^2+h(X,Y,Z).^2)+eps;

figure
scatter3(X(:),Y(:),Z(:),10*V(:),V(:))