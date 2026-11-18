// Gmsh project created on Wed Oct 01 18:31:38 2025
SetFactory("OpenCASCADE");

shell = 32;
d = shell/6;
x0 = 4.6;
y0 = 54.4/2;

ShapeFromFile("Ship.step")

Sphere(2) = {x0, y0, 0, shell, -Pi/2, Pi/2, 2*Pi};
Sphere(3) = {x0, y0, 0, shell+d, -Pi/2, Pi/2, 2*Pi};
Sphere(4) = {x0, y0, 0, shell+2*d, -Pi/2, Pi/2, 2*Pi};
Sphere(5) = {x0, y0, 0, shell+3*d, -Pi/2, Pi/2, 2*Pi};

Coherence;

c = newc;
p = newp;
Point(p) = {x0,y0-25,-3};
Point(p+1) = {x0,y0+25,-3};
Curve(c) = {p,p+1};
Line{c} In Volume{2};
Point(p+2) = {x0,y0+12,3};
Point(p+3) = {x0,y0+2,3};
Curve(c+1) = {p+2,p+3};
Line{c+1} In Volume{2};

Physical Volume("M:Hull") = 1;
Physical Volume("M:Air") = 2;
Physical Volume("M:Layer1") = 3;
Physical Volume("M:Layer2") = 4;
Physical Volume("M:Layer3") = 5;

Physical Surface("BD:Outer") = 60;

Characteristic Length{:} = d/2;
Characteristic Length{PointsOf{Volume{1};}} = 0.1;
Characteristic Length{PointsOf{Line{c};}} = 0.5;
Characteristic Length{PointsOf{Line{c+1};}} = 0.5;

Mesh 3;

Mesh.MshFileVersion = 2;

Save "ship.msh"; 