// Gmsh project created on Mon Oct  6 20:25:20 2025
SetFactory("OpenCASCADE");


shell = 35;
d = shell/6;
x0 = 4.6;
y0 = 54.4/2;
bd = 1;

ShapeFromFile("nutshell.step")

Sphere(2) = {0, 0, 0, shell, -Pi/2, Pi/2, 2*Pi};
Sphere(3) = {0, 0, 0, shell+d, -Pi/2, Pi/2, 2*Pi};
Sphere(4) = {0, 0, 0, shell+2*d, -Pi/2, Pi/2, 2*Pi};
Sphere(5) = {0, 0, 0, shell+3*d, -Pi/2, Pi/2, 2*Pi};

Coherence;

v = newv;
Box(v) = {8.5+bd, -25-bd , 2+bd, -17-2*bd, 50+2*bd, -12-2*bd};
Box(v+1) = {3.6-bd, -20+bd, -8+bd, -7.2+2*bd, 36.5-2*bd, 8-2*bd};
Delete{Volume{v,v+1};}

Surface{57:62} In Volume{2};
Surface{63:68} In Volume{3};

//c = newc;
//p = newp;
//Point(p) = {x0,y0-25,-3};
//Point(p+1) = {x0,y0+25,-3};
//Curve(c) = {p,p+1};
//Line{c} In Volume{2};
//Point(p+2) = {x0,y0+12,3};
//Point(p+3) = {x0,y0+2,3};
//Curve(c+1) = {p+2,p+3};
//Line{c+1} In Volume{2};

Physical Volume("M:Hull") = 1;
Physical Volume("M:Air") = {2,3};
Physical Volume("M:Layer1") = 4;
Physical Volume("M:Layer2") = 5;
Physical Volume("M:Layer3") = 6;

Physical Surface("BD:Outer") = 56;

Characteristic Length{:} = d/2;
Characteristic Length{PointsOf{Volume{1};}} = 0.25;
Characteristic Length{PointsOf{Surface{1:32};}} = 0.1;
Characteristic Length{PointsOf{Surface{57:62};}} = 1.5;
Characteristic Length{PointsOf{Surface{63:68};}} = 1.25;
//Characteristic Length{PointsOf{Line{c+1};}} = 0.5;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "nutshell.msh"; 