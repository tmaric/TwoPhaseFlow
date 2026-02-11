// Gmsh project created for piston radiation (axisymmetric wedge)
SetFactory("OpenCASCADE");
Geometry.OCCSewFaces = 0;

// Acoustic parameters
f = 6000;
c = 343;
lambda = c/f;

// Mesh resolution
N = 60; // cells per wavelength
gr=1.00;
gs = lambda/N;

// Geometry parameters
// piston radius
ap = 0.005;
// acoustic domain radius (PML starts here)
R0 = 0.2;
// PML thickness
L = 0.04;
// outer radius
R = R0 + L;

// Wedge angle (axisymmetric)
angle = 2*Pi/180;
rotateHalf = -angle/2;

// Points (x, y, z) with y = 0 plane
Point(1) = {0, 0, 0, gs};
Point(2) = {R, 0, 0, gs};
Point(3) = {0, R, 0, gs};

// Curves
Line(1) = {1, 2}; // piston (z=0, r<=ap)
Line(2) = {1, 3}; // symmetry axis (r=0), outer segment
Circle(3) = {2, 1, 3}; // outer spherical boundary (quarter circle)

// Surfaces: inner fan (unstructured) + outer quad ring (structured)
Curve Loop(1) = {1, 3, -2};
Plane Surface(1) = {1};

// Center the wedge around the symmetry plane
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


// Boundary patches (explicit IDs to avoid multi-physical surface merges)
Physical Surface("front") = {1};
Physical Surface("back") = {4};
Physical Surface("bottom") = {2};
Physical Surface("openAir") = {3};

Physical Volume("internal") = {1};

Mesh.ElementOrder = 1;
Mesh 3;
Save "pistonRadiation.msh";
