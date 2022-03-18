// Gmsh project created on Wed Mar 16 13:48:07 2022
SetFactory("OpenCASCADE");
General.NumThreads = 12;

Merge "aircoil-Body.brep";

Cylinder(2) = {0,0,-50,0,0,100,30,2*Pi};

Coherence;

Physical Volume("M:Coil") = {1,2,3};
Physical Volume("M:Air") = {4};

Physical Surface("BD:Upper") = 13;
Physical Surface("BD:Lower") = 11;

Characteristic Length{PointsOf{Volume{4};}} = 5;
Characteristic Length{PointsOf{Volume{1,2,3};}} = .5;

Mesh.Algorithm3D = 10; // For HXT

Mesh 3;

Mesh.MshFileVersion = 2;

Save "aircoil.msh";