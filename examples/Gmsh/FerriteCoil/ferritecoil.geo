// Gmsh project created on Thu Mar 17 13:17:39 2022
SetFactory("OpenCASCADE");

General.NumThreads = 12;

Merge "ferritecoil.brep";

Cylinder(6) = {0,0,-50,0,0,100,30,2*Pi};
Coherence;

Physical Volume("M:Coil") = {3:7};
Physical Volume("M:Core") = 2;
Physical Volume("M:Air") = 8;

Physical Surface("BD:UpperPort") = 22;
Physical Surface("BD:LowerPort") = 17;

Characteristic Length{PointsOf{Volume{8};}} = 5;
Characteristic Length{PointsOf{Volume{3:7};}} = .5;
Characteristic Length{PointsOf{Volume{2};}} = 1;

Mesh.Algorithm3D = 10;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "ferritecoil.msh";
