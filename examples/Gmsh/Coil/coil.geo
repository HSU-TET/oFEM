// Gmsh project created on Tue Mar 15 14:52:45 2022
SetFactory("OpenCASCADE");
General.NumThreads = 10;
Merge "Coil.brep";//+
Cylinder(3) = {0, 0, 0, 0, 0, 100, 25, 2*Pi};
Coherence;
//Recursive Delete{Volume{4,5};}

Physical Surface("BD:Upper") = 45;
Physical Surface("BD:Lower") = 43;
Physical Surface("BD:PEC") = {21,22,47,48};

Physical Volume("M:Ferrite") = 2;
Physical Volume("M:Coil") = {3,4,5};
Physical Volume("M:Air") = 6;

Characteristic Length{PointsOf{Volume{6};}} = 5;
Characteristic Length{PointsOf{Volume{3,4,5};}} = .5;
Characteristic Length{PointsOf{Volume{2};}} = 1;

Mesh.Algorithm3D = 10; // For HXT

Mesh 3;

Mesh.MshFileVersion = 2;

Save "Coil.msh";