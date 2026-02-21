// Gmsh project created for piston radiation (axisymmetric wedge)
SetFactory("OpenCASCADE");
Geometry.OCCSewFaces = 0;

f = 6000;                    // driven frequency
c = 343;                         // air sound speed
lambda = c/f;                       // wave length

N = 50;      // cells per wavelength
gs = lambda/N;                      // cell size
gr = 1.0;                   // transfinite progression ratio

R0 = 0.2;                  // acoustic domain radius (PML starts here)
L = 0.08;                      // PML thickness
R = R0 + L;                         // outer radius

// Wedge angle (axisymmetric)
angle = 2*Pi/180;
rotateHalf = -angle/2;

Point(1) = {0, 0, 0, gs};
Point(2) = {R, 0, 0, gs};
Point(3) = {0, R, 0, gs};

Line(1) = {1, 2};
Line(2) = {1, 3};
Circle(3) = {2, 1, 3};

Curve Loop(1) = {1, 3, -2};
Plane Surface(1) = {1};

Rotate{ {0,1,0}, {0,0,0}, rotateHalf } { Surface{1}; }

Transfinite Curve {1,3} = Round(R/gs) Using Progression gr;
Transfinite Curve {2} = Round(0.4*R/gs);
Transfinite Surface {1};
Recombine Surface {1};

out1[] = Extrude {{0,1,0}, {0,0,0}, angle} {
  Surface{1};
  Layers{1};
  Recombine;
};

Physical Surface("front") = {1};
Physical Surface("back") = {4};
Physical Surface("bottom") = {2};
Physical Surface("openAir") = {3};

Physical Volume("internal") = {1};

Mesh.ElementOrder = 1;
Mesh 3;
Save "pistonRadiation.msh";
