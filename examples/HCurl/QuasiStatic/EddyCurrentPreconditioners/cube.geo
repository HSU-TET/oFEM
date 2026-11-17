// Cube [-pi/2,pi/2]^3 with PEC boundary (copy of RectangularCavityConvergence/rect.geo)
SetFactory("OpenCASCADE");

DefineConstant[ h = 2 ];   // refinement (mesh size pi/(2h))

x = 1;
y = 1;
z = 1;

Box(1) = {-x/2*Pi,-y/2*Pi,-z/2*Pi,x*Pi,y*Pi,z*Pi};

Physical Volume("M:Domain") = 1;
Physical Surface("BD:PEC") = {1:6};

Characteristic Length{PointsOf{Volume{1};}} = Pi/(2*h);

//Mesh 3;

Mesh.MshFileVersion = 2;

//Save "rect.msh";