// Channel with a cylinder (geometry of the DFG flow around a cylinder benchmark)
// Mesh with:  gmsh channel.geo -setnumber h 1 -2 -o channel.msh
// h refines the mesh (mesh size ~ 1/h), xc and yc move the cylinder.
SetFactory("OpenCASCADE");

DefineConstant[ h = 1, xc = 0.2, yc = 0.2 ];   // refinement, cylinder center

L  = 2.2;   // channel length
H  = 0.41;  // channel height
r  = 0.05;  // cylinder radius

Rectangle(1) = {0, 0, 0, L, H};
Disk(2) = {xc, yc, 0, r};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

e = 1e-6;
inlet()  = Curve In BoundingBox{-e, -e, -e, e, H+e, e};
outlet() = Curve In BoundingBox{L-e, -e, -e, L+e, H+e, e};
bottom() = Curve In BoundingBox{-e, -e, -e, L+e, e, e};
top()    = Curve In BoundingBox{-e, H-e, -e, L+e, H+e, e};
cyl()    = Curve In BoundingBox{xc-r-e, yc-r-e, -e, xc+r+e, yc+r+e, e};

Physical Surface("M:Domain") = {3};
Physical Curve("BD:Inlet") = {inlet()};
Physical Curve("BD:Outlet") = {outlet()};
Physical Curve("BD:Wall") = {bottom(), top()};
Physical Curve("BD:Cylinder") = {cyl()};

// mesh size: lc in the channel, lc/5 at the cylinder
lc = 0.04/h;
Field[1] = Distance;
Field[1].CurvesList = {cyl()};
Field[1].Sampling = 200;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc/5;
Field[2].SizeMax = lc;
Field[2].DistMin = r/4;
Field[2].DistMax = 6*r;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.MshFileVersion = 2;
