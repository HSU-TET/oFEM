// Gmsh project created on Mon Oct 21 11:21:59 2019
SetFactory("OpenCASCADE");
Coherence; 

width =0.02;
hight = 0.01; 
depth = 0.01; 

Box(1) = {0,0,0, width/2, hight, depth}; 
Box(2) = {width/2, 0, 0, width/2, hight, depth }; 
Coherence;
Physical Volume("M:LeftCap") = 1; 
Physical Volume("M:RightCap") = 2; 
Coherence;
Physical Surface("BD:LeftDirichlet")= 1; 
Physical Surface("BD:RightDirichlet")= 7; 
Coherence;

Characteristic Lenght{:} = hight/10;


Mesh 3; 

Mesh.MshFileVersion = 2; 
Save "planarCapacitor3D.msh";

