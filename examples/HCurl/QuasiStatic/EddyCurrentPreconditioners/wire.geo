// Straight copper wire in an air box, PEC on the outer boundary.
// Mesh with:  gmsh wire.geo -setnumber h 1 -3 -o wire.msh
SetFactory("OpenCASCADE");

DefineConstant[ h = 1 ];   // refinement (mesh size ~ 1/h)

W = 0.1;     // box width
L = 0.1;     // box (and wire) length
r = 0.001;   // wire radius

Box(1) = {-W/2, -W/2, 0, W, W, L};
Cylinder(2) = {0, 0, 0, 0, 0, L, r};
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

e = 1e-6;
wire() = Volume In BoundingBox{-r-e, -r-e, -e, r+e, r+e, L+e};
all() = Volume{:};
air() = all();
air() -= wire();

Physical Volume("M:Wire") = {wire()};
Physical Volume("M:Air") = {air()};
Physical Surface("BD:PEC") = CombinedBoundary{ Volume{all()}; };

// mesh size: lc in the air, r/2 at the wire
lc = W/(6*h);
Field[1] = Distance;
Field[1].SurfacesList = {Boundary{ Volume{wire()}; }};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = r/(2*h);
Field[2].SizeMax = lc;
Field[2].DistMin = r;
Field[2].DistMax = W/4;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.MshFileVersion = 2;
