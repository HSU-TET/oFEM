// Gmsh project created on Thu Dec 12 17:24:03 2019
SetFactory("OpenCASCADE");

r_skin=0.06 ;
r_bone=0.008 ; 
r_muscle=0.04;
r_fat=0.055; 
w=0.004;
h=0.001;

Disk(1)= {0,0,0,r_skin};
Disk(2)={0,0,0,r_fat};
//Disk(3)={0,0,0,r_muscle};
//Disk(4)={0.01,0,0,r_bone}; 
//Disk(5)={-0.01,0,0,r_bone};
Rectangle(100)= {0.058,0.01,0,w,h};
Rectangle(200)= {0.058,-0.01,0,w,h};
Coherence; 

Physical Surface("M:skin")= {3,4,5}; 
//Physical Surface("M:muscle")= {10}; 
Physical Surface("M:fat")= {2}; 
//Physical Surface("M:bone")= {4,5}; 

Physical Line("BD:leftElectrode")={11};
Physical Line("BD:rightElectrode")={12};

//Delete{Line{2,3,4,6,7,8};}
Delete{Surface{6,7};}
p = newp;
Point(p) = {0,0,0};
Point{p} In Surface{2};

Characteristic Length{PointsOf{Surface{2,3};}} = 0.0005;
Characteristic Length{PointsOf{Surface{1};} = 0.01;
Characteristic Length{p} = 0.001;

//Characteristic Length {:} = 0.001;
Mesh 2; 
Mesh.MshFileVersion = 2; 

Save "Simple2D.msh";
