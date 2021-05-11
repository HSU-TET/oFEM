// Gmsh project created on Sat Dec 07 12:46:20 2019
SetFactory("OpenCASCADE");
Coherence; 


h = 0.33; 
thickness_skin= 0.004; 
thickness_fat=0.01; 

r_skinmax = 0.6/(2*3.1416);
r_skinmedium = 0.52/(2*3.1416);
r_skinmin = 0.4/(2*3.1416);

r_fatmax = r_skinmax - thickness_skin;
r_fatmedium = r_skinmedium - thickness_skin;
r_fatmin = r_skinmin - thickness_skin;

r_musclemax = r_fatmax - thickness_fat; 
r_musclemedium = r_fatmedium - thickness_fat;
r_musclemin = r_fatmin - thickness_fat; 

r_bone = 0.02;
r_vein = 0.005; 

// Muscle
Circle(13)= {0,0,0,r_musclemin};
Circle(14)={0,0,0.13,r_musclemedium};
Circle(15)={0,0,0.33,r_musclemax};
Curve Loop(16) = {13};
Curve Loop(17) = {14};
Curve Loop(18) = {15};
Ruled ThruSections(3) = {16,17,18}; 


// Fat
Circle(7)= {0,0,0,r_fatmin};
Circle(8)={0,0,0.13,r_fatmedium};
Circle(9)={0,0,0.33,r_fatmax};
Curve Loop(10) = {7};
Curve Loop(11) = {8};
Curve Loop(12) = {9};
Ruled ThruSections(2) = {10,11,12}; 
Coherence; 

// Skin
Circle(1)= {0,0,0,r_skinmin};
Circle(2)={0,0,0.13,r_skinmedium};
Circle(3)={0,0,0.33,r_skinmax};
Curve Loop(4) = {1};
Curve Loop(5) = {2};
Curve Loop(6) = {3};
Ruled ThruSections(1) = {4,5,6}; 
Coherence; 


//Bone
Cylinder(6) = {0           ,0    ,0,0,0,h,r_bone};

//Blood
Cylinder(7) = {r_bone+0.005,0    ,0,0,0,h,r_vein};
Cylinder(8) = {r_bone+0.02 ,0.005,0,0,0,h,r_vein};
Coherence;  

//Electodes. 
Point(20) ={0.08355688845,0.05/2,0.2}; 
Point(21) ={0.08355688845,-0.05/2,0.2}; 
diff= 0.08355688845-0.081;
Cylinder(20) ={0.081, 0.05/2,0.2, diff + 0.001, 0, 0, 0.003 }; 
Coherence; 
Cylinder(21) ={0.081,-0.05/2,0.2, diff + 0.001, 0, 0, 0.003 }; 
Coherence; 


Physical Surface("BD:leftElectrode")= 34; 
Physical Surface("BD:rightElectrode") =43; 
Physical Surface("BD:UpperCircle")={8,25,20,26,40,14}; 
Physical Surface("BD:LowerCircle")={22,21,15,6,38}; 

Physical Volume("M:Skin") = {13,11,14}; 
Physical Volume("M:Fat") = 4; 
Physical Volume("M:Muscle") =9 ; 
Physical Volume("M:Blood") ={7,8};
Physical Volume("M:Bone") = 6; 
// Volume 12,15 Electode 

Characteristic Length {:} = 0.005;
Mesh 3; 
Mesh.MshFileVersion = 2; 

Save "leg3D.msh";