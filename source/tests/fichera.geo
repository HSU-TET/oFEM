// Gmsh project created on Thu Jun 11 12:26:24 2020
SetFactory("OpenCASCADE");

l = 2;
b = 2;
h = 2;
lh = 1;
bh = 1;
hh = 1;

Box(1) = {-l/2,-h/2,-b/2,l,h,b};
Box(2) = {0,0,0,lh,hh,bh};

BooleanDifference {Volume{1};Delete;}{Volume{2};Delete;}

Physical Volume("M:Domain") = 1;
Physical Surface("BD:Boundary") = {7,9,11:17};

Characteristic Length{:} = 0.2;
Characteristic Length{10} = 0.005;
Characteristic Length{9,11:14,16} = 0.05;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "fichera.msh";