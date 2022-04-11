// Gmsh project created on Thu Apr 07 13:51:29 2022
SetFactory("OpenCASCADE");

h=h;

Rectangle(1) = {-Pi/2,-Pi/2,0,Pi,Pi};

Physical Surface("M:Domain") = 1;
Physical Line("BD:PEC") = {1:4};

Characteristic Length{:} = Pi/(4*h);

Mesh 2;

Mesh.MshFileVersion = 2;

Save "rect.msh";