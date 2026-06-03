// Gmsh project created on Tue Apr 15 14:40:01 2025
SetFactory("OpenCASCADE");
// frequency
f=25250;
// sound speed
c=343;
// wavelength 
l=c/f;
// distance 
D = 0.00744;//1.586*l;
// mesh size
N=200;
gs=l/N;
lc=l/N;
lc1=0.4*l/N;

//z-coordinate of sphere center
s = 0.004;
//sphere radius
r = 0.001;
R = 2*r;
a=60*Pi/180;
b=70*Pi/180;

//geometry
rT = 0.01; // radius of transducer
rR = 0.019; // radius of reflector
bu = 0.01; // air buffer
hS = 0.01; // height of reflector/transducer stand

// points
Point(1) = {0, 0, 0, lc};
Point(2) = {rT, 0, 0, lc};
Point(3) = {rR, 0, 0, lc};
Point(4) = {rR, -hS,0, lc};
Point(5) = {rR+bu,  -hS, 0,lc};
Point(6) = {rR+bu,  0, 0,lc};
Point(7) = {rR+bu, D, 0,lc};
Point(8) = {rR+bu, D+hS,0, lc};
Point(9) = {rR, D+hS,0,lc};
Point(10) = {rT, D+hS,0, lc};
Point(11) = {rT, D,0,lc};
Point(12) = {0, D,0, lc};
Point(13) = {rR, D,0, lc};
Point(14) = {0, s-r,0, lc};
Point(15) = {0, s,0, lc};
Point(16) = {0, s+r,0, lc};
Point(17) = {R*Sin(a), D, 0, lc};
Point(18) = {R*Sin(a), 0, 0, lc};
Point(19) = {0, s+0.85*R, 0, lc};
Point(20) = {0, s-0.85*R, 0, lc};
Point(21) = {r*Sin(b), s+r*Cos(b), 0, lc};
Point(22) = {r*Sin(b), s-r*Cos(b), 0, lc};
Point(23) = {rT, s+R*Cos(a), 0, lc};
Point(24) = {rT, s-R*Cos(a), 0, lc};
Point(25) = {R*Sin(a), s+R*Cos(a), 0, lc};
Point(26) = {R*Sin(a), s-R*Cos(a), 0, lc};
Point(27) = {rR, s-R*Cos(a), 0, lc};
Point(28) = {rR, s+R*Cos(a), 0, lc};
Point(29) = {rR+bu, s-R*Cos(a), 0, lc};
Point(30) = {rR+bu, s+R*Cos(a), 0, lc};


Line(1) = {1,18};
Line(2) = {18,2};
Line(3) = {2,24};
Line(4) = {23,24};
Line(5) = {23,11};
Line(6) = {11,17};
Line(7) = {17,12};
Line(8) = {12,19};
Line(9) = {19,16};
Circle(10) = {16,15,21};
Circle(11) = {21,15,22};
Circle(12) = {22,15,14};
Line(13) = {14,20};
Line(14) = {20,1};
Line(15) = {25,23};
Line(16) = {25,17};
Line(17) = {26,18};
Line(18) = {24,26};
Circle(19) = {25,15,26};
Line(20) = {21,25};
Line(21) = {22,26};
Line(22)= {19,25};
Line(23)= {20,26};
Line(24)= {3,4};
Line(25)= {4,5};
Line(26)= {5,6};

Line(28)= {7,8};
Line(29)= {8,9};
Line(30)= {9,10};
Line(31)= {10,11};
Line(32)= {2,3};

Line(34)= {13,11};
Line(35)= {3,6};
Line(36)= {7,13};
Line(37)= {13,9};

Line(38)= {3,27};
Line(39)= {27,28};
Line(40)= {28,13};
Line(41)= {28,23};
Line(42)= {27,24};

Line(43)= {6,29};
Line(44)= {29,30};
Line(45)= {30,7};
Line(46)= {29,27};
Line(47)= {30,28};


Transfinite Line {25,35,36,29,46,47} = Round(bu/gs);
Transfinite Line {31,37,28,24,26} = Round(hS/gs);
Transfinite Line {30,32,34,42,41} = Round((rR-rT)/gs);
Transfinite Line {33,27} = Round(D/gs)-2;

Curve Loop(300)={24,25,26,-35};
Plane Surface(400)={300};

Curve Loop(302)={-36,28,29,-37};
Plane Surface(402)={302};

Curve Loop(304)={-34,37,30,31};
Plane Surface(404)={304};

Curve Loop(308)={-41,40,34,5};
Plane Surface(408)={308};

Curve Loop(306)={-47,45,36,-40};
Plane Surface(406)={306};

Curve Loop(307)={-42,39,41,4};
Plane Surface(407)={307};

Curve Loop(312)={-46,44,47,-39};
Plane Surface(412)={312};

Curve Loop(314)={32,38,42,-3};
Plane Surface(414)={314};

Curve Loop(316)={35,43,46,-38};
Plane Surface(416)={316};



Curve Loop(100)={1,-17,-23,14};
Plane Surface(200)={100};

Curve Loop(101)={17,2,3,18};
Plane Surface(201)={101};

Curve Loop(102)={19,-18,-4,-15};
Plane Surface(202)={102};

Curve Loop(103)={15,5,6,-16};
Plane Surface(203)={103};

Curve Loop(104)={22,16,7,8};
Plane Surface(204)={104};

Curve Loop(105)={10,20,-22,9};
Plane Surface(205)={105};

Curve Loop(106)={-20,11,21,-19};
Plane Surface(206)={106};

Curve Loop(108)={12,13,23,-21};
Plane Surface(208)={108};

Rotate{ {0,1,0}, {0,0,0}, - Pi*1/180} {Surface{:};}

Transfinite Line {52,54,57,59} = Round((rT-R*Sin(a))/gs);
Transfinite Line {63,60,58,82,84} = Round((D-(s+R*Cos(a)))/gs);
Transfinite Line {51,49,53,91,92} = Round((s-R*Cos(a))/gs);
Transfinite Curve {55,56,86,89,67} = Round(2*R*Cos(a)/gs);
Transfinite Curve {48,50,69,64,61,62} = 2*Round(R*Cos(a)/gs);
Transfinite Line {66,65,68,70} = Round((R-r)/gs)  Using Progression 1.1;

Transfinite Surface{:};
Recombine Surface{:};

surfaceVector[] = Extrude {{0,1,0}, {0,0,0}, Pi*2/180} {
	 Surface{:};
	 Layers{1};
	 Recombine;
	};


Physical Surface("transducer1") = {430, 434};

Physical Surface("transducer2") = {454};

Physical Surface("openAir") = {445,446,473,468,451,452,456,460};

Physical Surface("reflector1") = {421,470,417};

Physical Surface("reflector2") = {444};

Physical Surface("front") = {204,205,206,208,200,201,202,203,408,407,414,406,412,416,400,402,404};

Physical Surface("back") = {420,443,441,438,435,432,428,424,472,465,466,474,469,448,453,457,461};

Physical Surface("dropWall") = {442,439,436};

Physical Volume("internal")={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};


Mesh.ElementOrder = 1;
Mesh 3;
Save "levitatorSphereWedgeHex_refine.msh";




