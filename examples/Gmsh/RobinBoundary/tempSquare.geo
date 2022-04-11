// Gmsh project created on Thu Apr 01 11:52:23 2021
SetFactory("OpenCASCADE");

h = 10;

Rectangle(1) = {-h/2,-h/2,0,h,h};

Physical Surface("M:Domain") = 1;

Physical Line("BD:Left") = 4;
Physical Line("BD:Right") = 2;

Characteristic Length{:} = 0.5;

Mesh.MshFileVersion = 2;

Mesh 2;

Save "tempSquare.msh";