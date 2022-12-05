// Gmsh project created on Fri Nov 22 18:14:41 2019
SetFactory("OpenCASCADE");

Coherence;

width =0.04 ; 
height = 0.04; 

Rectangle(1) ={-width/2,-height/2,0, width/2, height}; 
Rectangle(2) ={0,-height/2,0, width/2, height};

Coherence;

Physical Line("BD:LeftPlate") = {4};
Physical Line("BD:RightPlate")= {6};

Physical Surface("M:Paper") = {1};
Physical Surface("M:Glas") = {2};

Characteristic Length{:}= height/h;

Mesh 2; 

Mesh.MshFileVersion = 2; 

Save "planarCapacitor.msh";

Geometry.LineNumbers = 1;
Geometry.Surfaces = 1;
Geometry.SurfaceNumbers = 1;