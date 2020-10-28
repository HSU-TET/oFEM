// Gmsh project created on Fri Nov 22 19:30:11 2019
SetFactory("OpenCASCADE");

Coherence;

length =1 ; 
higth = 0.1; 

Rectangle(1) ={0,0,0, length, higth}; 
//Rectangle(2) ={0,higth,0, length, higth};

Physical Line("BD:HighTemp") = {4};
Physical Line("BD:LowTemp")= {2};
Physical Line("BD:Boundary") = {1,3};
//Physical Line("BD:TextileBoundary")= {7};

Physical Surface("M:Steel") = {1};
//Physical Surface("M:Textile") = {2};

Characteristic Length{:}= 0.02;

Mesh 2; 

Mesh.MshFileVersion = 2; 

Save "HeatExample2d.msh";
