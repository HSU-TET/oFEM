// Gmsh project created on Fri Nov 22 18:14:41 2019
SetFactory("OpenCASCADE");

Coherence;

width =0.04 ; 
higth = 0.04; 

Rectangle(1) ={0,0,0, width/2, higth}; 
Rectangle(2) ={width/2,0,0, width/2, higth};

Coherence;

Physical Line("BD:left") = {4};
Physical Line("BD:right")= {6};

Physical Surface("M:Paper") = {1};
Physical Surface("M:Glas") = {2};

Characteristic Length{:}= 0.01;

Mesh 2; 

Mesh.MshFileVersion = 2; 

Save "CapacitorExample2d.msh";



