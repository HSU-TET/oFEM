// Gmsh project created on Mon Jan 27 14:05:21 2025
SetFactory("OpenCASCADE");

shell = 0.25;
d = shell/6;

Torus(1) = {0, 0, 0, 0.1, 0.05, 2*Pi};

Sphere(2) = {-0, -0, 0, shell, -Pi/2, Pi/2, 2*Pi};
Sphere(3) = {-0, -0, 0, shell+d, -Pi/2, Pi/2, 2*Pi};
Sphere(4) = {-0, -0, 0, shell+2*d, -Pi/2, Pi/2, 2*Pi};
Sphere(5) = {-0, -0, 0, shell+3*d, -Pi/2, Pi/2, 2*Pi};

Coherence;
Box(6) = {-0.125,-0.125,0.1,0.25,0.25,0.04};
Box(7) = {-0.125,-0.125,-0.1,0.25,0.25,-0.04};
Torus(8) = {0,0,0,0.1,0.06,2*Pi};
Torus(9) = {0,0,0,0.1,0.025,2*Pi};

Delete{Volume{6,7,8,9};}
Surface{6:17} In Volume{2};
Surface{18} In Volume{2};
Surface{19,20} In Volume{1};

Physical Volume("M:Coil") = 1;
Physical Volume("M:Air") = 2;
Physical Volume("M:Layer1") = 3;
Physical Volume("M:Layer2") = 4;
Physical Volume("M:Layer3") = 5;

Physical Surface("BD:Outer") = 5;

Physical Surface("ROI:Upper") = 10;
Physical Surface("ROI:Lower") = 17;

Characteristic Length{:} = d/2;
Characteristic Length{PointsOf{Volume{1};}} = d/16;
Characteristic Length{PointsOf{Surface{6:17};}} = 0.02;
Characteristic Length{PointsOf{Surface{18};}} = d/8;
Characteristic Length{PointsOf{Surface{19};}} = d/8;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "coilInAir.msh"; 