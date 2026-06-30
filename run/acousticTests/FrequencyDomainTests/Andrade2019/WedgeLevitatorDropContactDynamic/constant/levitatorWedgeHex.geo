// Gmsh project created for Andrade2019 levitator wedge case
SetFactory("OpenCASCADE");

// frequency
f = 25250;
// sound speed
c = 343;
// wavelength
l = c/f;

// distance
H = 0.00764;
// mesh size controls
lc = 1.0;
N = 200;
gs = l/N;

// geometry
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
Point(7) = {rR+bu, H, 0,lc};
Point(8) = {rR+bu, H+hS,0, lc};
Point(9) = {rR, H+hS,0,lc};
Point(10) = {rT, H+hS,0, lc};
Point(11) = {rT, H,0,lc};
Point(12) = {0, H,0, lc};
Point(13) = {rR, H,0, lc};

// lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,1};
Line(13) = {2,11};
Line(14) = {3,13};
Line(15) = {13,11};
Line(16) = {13,9};
Line(17) = {13,7};
Line(18) = {6,3};

Curve Loop(19)={1,13,11,12};
Plane Surface(20)={19};

Curve Loop(21)={2,14,15,-13};
Plane Surface(22)={21};

Curve Loop(23)={3,4,5,18};
Plane Surface(24)={23};

Curve Loop(25)={-18,6,-17,-14};
Plane Surface(26)={25};

Curve Loop(27)={17,7,8,-16};
Plane Surface(28)={27};

Curve Loop(29)={16,9,10,-15};
Plane Surface(30)={29};

Rotate{ {0,1,0}, {0,0,0}, - Pi*0.5/180} {Surface{20, 22, 24, 26 , 28 , 30};}

Transfinite Line {19, 21} = Round(rT/gs);
Transfinite Line {23, 25, 35} = Round((rR-rT)/gs);
Transfinite Line {22,20,24,30}=Round(H/gs);
Transfinite Line {27,29,31,33}=Round(bu/gs);
Transfinite Line {32,34,36,26,28}=Round(hS/gs);

Transfinite Surface "*";
Recombine Surface "*";

surfaceVector[] = Extrude {{0,1,0}, {0,0,0}, Pi*1.0/180} {
     Surface{20, 22, 24, 26 , 28 , 30};
     Layers{1};
     Recombine;
    };

Physical Surface("transducer1") = {33};
Physical Surface("transducer2") = {52};
Physical Surface("openAir") = {51, 48,47,44,41, 40};
Physical Surface("reflector1") = {35,31};
Physical Surface("reflector2") = {39};
Physical Surface("front") = {20, 22, 24, 26, 28, 30};
Physical Surface("back") = {34, 38, 53, 50, 46, 43};

Physical Volume("internal")={1,2,3,4,5, 6 };

Mesh.ElementOrder = 1;
Mesh 3;
Save "levitatorWedgeHex.msh";
