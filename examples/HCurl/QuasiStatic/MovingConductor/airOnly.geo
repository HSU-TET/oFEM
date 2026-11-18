// Gmsh project created on Mon Jan 27 14:05:21 2025
SetFactory("OpenCASCADE");

shell = 0.25;
d = shell/6;

Sphere(1) = {-0, -0, 0, shell, -Pi/2, Pi/2, 2*Pi};
Sphere(2) = {-0, -0, 0, shell+d, -Pi/2, Pi/2, 2*Pi};
Sphere(3) = {-0, -0, 0, shell+2*d, -Pi/2, Pi/2, 2*Pi};
Sphere(4) = {-0, -0, 0, shell+3*d, -Pi/2, Pi/2, 2*Pi};

Coherence;

Physical Volume("M:Air") = 1;
Physical Volume("M:Layer1") = 2;
Physical Volume("M:Layer2") = 3;
Physical Volume("M:Layer3") = 4;

Physical Surface("BD:Outer") = 4;

Characteristic Length{:} = d;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "airOnly.msh"; 