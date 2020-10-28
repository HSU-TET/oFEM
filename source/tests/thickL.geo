// Gmsh project created on Thu Jun 11 13:40:14 2020
SetFactory("OpenCASCADE");

l = 2;
b = 2;
lh = 1;
bh = 1;

h = 1;

Box(1) = {-l/2,0,-b/2,l,h,b};
Box(2) = {0,0,0,lh,h,bh};

BooleanDifference {Volume{1};Delete;}{Volume{2};Delete;}

Physical Volume("M:Domain") = 1;
Physical Surface("BD:Boundary") = {7,11:17};

Characteristic Length{:} = 0.25;
Characteristic Length{PointsOf{Curve{16};}} = 0.05;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "thickL.msh";