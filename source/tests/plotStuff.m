%f = @(x,y,z) [1-2.*x-y-z,-x,-x];
%f = @(u,v,w) [-v-u.*v+v.^2+v.*w,-u+u.^2-u.*v+u.*w,-2*u.*v];
f = @(u,v,w) [3*u.^2.*v-6*u.*v.^2+v.^3,u.^3-6*u.^3.*v+3*u.*v.^2];

[X,Y,Z] = meshgrid(linspace(0,1,11));

X(X+Y+Z>1) = NaN;
Y(X+Y+Z>1) = NaN;
Z(X+Y+Z>1) = NaN;

% X = X(:,:,1);
% Y = Y(:,:,1);
% Z = Z(:,:,1);

V = f(X,Y,Z);

V1 = V(:,1:size(X,1),:);
V2 = V(:,1+size(X,1):2*size(X,1),:);
%V3 = V(:,1+2*size(X,1):3*size(X,1),:);
V3 = zeros(size(V1));

quiver3(X,Y,Z,V1,V2,V3)