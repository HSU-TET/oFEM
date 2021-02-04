// Gmsh project created on Mon Jun 22 17:53:52 2020
SetFactory("OpenCASCADE");

phi = Pi/4;
rc = 1e-3;
rA = 5e-3;
l = 5e-3;

Cylinder(1) = {0,0,0,0,0,l,rc,phi};
Cylinder(2) = {0,0,0,0,0,l,rA,phi};

Coherence;

Physical Volume("M:Conductor") = 1;
Physical Volume("M:Air") = 2;

Physical Surface("BD:PMC") = {3,5,9,10};
Physical Surface("BD:PEC") = {1,2,4,7,8};

Characteristic Length{:} = 5e-4;
Characteristic Length{PointsOf{Volume{1};}} = 2e-4;
Characteristic Length{PointsOf{Surface{6};}} = 1e-4;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "pie.msh";