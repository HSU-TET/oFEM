// Gmsh project created on Thu Jun 11 15:22:55 2020
SetFactory("OpenCASCADE");

r = 0.001;
l = 0.05;
phi = 2*Pi;
rAir = 0.01;

Cylinder(1) = {0,0,0,0,0,l,r};
Cylinder(2) = {0,0,0,0,0,l,rAir};

Coherence;

Physical Volume("M:Conductor") = 1;
Physical Volume("M:Air") = 2;

Physical Surface("BD:PEC") = {2,3,5,6};

Characteristic Length{PointsOf{Volume{2};}} = 0.001;
Characteristic Length{PointsOf{Volume{1};}} = 0.0005;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "testWire.msh";