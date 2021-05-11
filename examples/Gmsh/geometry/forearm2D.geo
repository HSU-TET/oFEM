// Gmsh project created on Thu Dec 12 17:24:03 2019
SetFactory("OpenCASCADE");

r_skin=0.06 ;
r_bone=0.008 ; 
r_muscle=0.04;
r_fat=0.05; 
w=0.004;
h=0.001;

Disk(1)= {0,0,0,r_skin};
Disk(2)={0,0,0,r_fat};
Disk(3)={0,0,0,r_muscle};
Disk(4)={0.01,0,0,r_bone}; 
Disk(5)={-0.01,0,0,r_bone};
Rectangle(100)= {0.058,0.01,0,w,h};
Rectangle(200)= {0.058,-0.01,0,w,h};
Coherence; 

Physical Surface("M:skin")= {6,7,8}; 
Physical Surface("M:muscle")= {10}; 
Physical Surface("M:fat")= {9}; 
Physical Surface("M:bone")= {4,5}; 

Physical Line("BD:leftElectrode")= {12};
Physical Line("BD:rightElectrode")= {11};

Delete{Line{2,3,4,6,7,8};}
Delete{Surface{11,12};}

Characteristic Length {:} = 0.001;
Mesh 2; 
Mesh.MshFileVersion = 2; 

Save "forearm2D.msh";
