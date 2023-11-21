// Define bounding cylinder to set up poly mesh for rising bubble
// test case.
// The axis of the cylinder points in y-direction
pi = 3.14159265359;

// Denner Setup
// Radius and height of cylinder
//r = 0.05;
//h = 0.14;

// old lent setup
//rDroplet = 0.5*6.94975e-03;
//r = 1.5e-02;
//h = 5e-02;
//es= 0.001;

rDroplet = 0.2;
r = 0.5;
h = 2;
es = 0.05;

// Define bottom centre of cylinder
x = r;
y = 0;
z = r;

// Define element size


// Bottom centre point
Point(1) = {x, y, z};

// Support points in x-z plane with y=0
Point(2) = {x, y, z-r};
Point(3) = {x+r, y, z};
Point(4) = {x, y, z+r};
Point(5) = {x-r, y, z};

// Support points in x-z plane with y=h
Point(6) = {x, y+h, z-r};
Point(7) = {x+r, y+h, z};
Point(8) = {x, y+h, z+r};
Point(9) = {x-r, y+h, z};

// Define bottom circle arcs
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Define top circle arcs
Point(10) = {x, y+h, z};

Circle(5) = {9, 10, 8};
Circle(6) = {8, 10, 7};
Circle(7) = {7, 10, 6};
Circle(8) = {6, 10 ,9};

// Define cylinder bottom segment
Line Loop(1) = {1, 2, 3, 4};

// Define cylinder top segment
Line Loop(2) = {5, 6, 7, 8};

// Define halfs of nappe
Line(9) = {3, 7};
Line(10) = {4, 8};
Line(11) = {5, 9};
Line(12) = {2, 6};

Line Loop(3) = {12, -7, -9, -1};
Line Loop(4) = {9, -6, -10, -2};
Line Loop(5) = {10, -5, -11, -3};
Line Loop(6) = {11, -8, -12, -4};

// Define cylinder surface patches
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};

// Transformation to obtain different axis orientation. Here, the cylinder
// is rotated around the x-Axis by pi/2, so that its axis points in z-direction
Rotate {{1, 0, 0}, {0, 0, 0}, 0.5*pi}{Surface{1,2,3,4,5,6};}
Translate {0, 2*r, 0}{Surface{1,2,3,4,5,6};}

// Combine to cylinder surface
// Only mesh one patch at a time with gmsh if you want do keep different patches
// with cfMesh
//Physical Surface(1) = {1, 2, 3, 4, 5, 6};
// Bottom patch
//Physical Surface(1) = {1};
// Top patch
//Physical Surface(2) = {2};
// Mantle patch
Physical Surface(3) = {3, 4, 5, 6};

Mesh.CharacteristicLengthMax = es;
Mesh.RemeshAlgorithm = 1; // Choose automatically
Mesh.RemeshParametrization = 0; // harmonic
Mesh.Algorithm = 6; // Frontal
