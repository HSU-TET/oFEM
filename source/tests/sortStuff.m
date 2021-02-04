% Generate a simple Test Mesh, to check geometry functionality
close all
% Instantiate geometry
mesh = ofem_v2.Geometry;

% Fill basic data
mesh.dim = 3;
mesh.type = 'tet';

%% Use a [0,1]^3 Cube, here with 6 tets
%x = 0:0.5:1;
%y = 0:0.5:1;
%z = 0:0.5:1;

%% Switch to cavity Jin page 272
% x = linspace(0,0.01,11);
% y = linspace(0,0.005,7);
% z = linspace(0,0.0075,11);

disc = 5;

x = linspace(0,pi,disc);
y = linspace(0,pi,disc);
z = linspace(0,pi,disc);

[X,Y,Z] = meshgrid(x,y,z);

X = X(:);
Y = Y(:);
Z = Z(:);

TR = delaunayTriangulation(X,Y,Z);

co = ofem_v2.tools.matrixarray(reshape([X,Y,Z]',3,1,[]));

mesh.el = TR.ConnectivityList;
mesh.co = ofem_v2.tools.matrixarray(reshape(TR.Points',3,1,[]));

mesh.Nco = size(mesh.co,3);

mesh.Nint = size(mesh.el,1);


















