// Gmsh project created on Tue Mar 22 14:07:55 2022
SetFactory("OpenCASCADE");

x = 1;
y = 1;
z = 1;

Box(1) = {-x/2*Pi,-y/2*Pi,-z/2*Pi,x*Pi,y*Pi,z*Pi};

Physical Volume("M:Domain") = 1;
Physical Surface("BD:PEC") = {1:6};

Characteristic Length{PointsOf{Volume{1};}} = Pi/2;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "rect.msh";