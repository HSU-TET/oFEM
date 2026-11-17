// Gmsh project created on Thu Apr 01 12:51:27 2021
SetFactory("OpenCASCADE");

h = 10;

Box(1) = {-h/2,-h/2,-h/2,h,h,h};

Physical Volume("M:Domain") = 1;

Physical Surface("BD:Left") = 1;
Physical Surface("BD:Right") = 2;

Characteristic Length{:} = 0.5;

Mesh.MshFileVersion = 2;

Mesh 3;

Save "tempCube.msh";