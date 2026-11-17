// Two layers (paper below, textile above) between a low and a high
// concentration on the left and right ends.
// Mesh with:  gmsh Diffusion2d.geo -2 -o Diffusion2d.msh
SetFactory("OpenCASCADE");

L = 1;     // length
H = 0.1;   // height of each layer

Rectangle(1) = {0, 0, 0, L, H};   // paper
Rectangle(2) = {0, H, 0, L, H};   // textile
Coherence;

e = 1e-6;
left()   = Curve In BoundingBox{-e, -e, -e, e, 2*H+e, e};
right()  = Curve In BoundingBox{L-e, -e, -e, L+e, 2*H+e, e};
bottom() = Curve In BoundingBox{-e, -e, -e, L+e, e, e};
top()    = Curve In BoundingBox{-e, 2*H-e, -e, L+e, 2*H+e, e};

// Dirichlet boundaries first: where boundaries share nodes, the first
// one keeps them
Physical Curve("BD:LowCon") = {left()};
Physical Curve("BD:HighCon") = {right()};
Physical Curve("BD:PaperBoundary") = {bottom()};
Physical Curve("BD:TextileBoundary") = {top()};

Physical Surface("M:Paper") = {1};
Physical Surface("M:Textile") = {2};

Characteristic Length{:} = 0.02;

Mesh.MshFileVersion = 2;
