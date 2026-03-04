// 3D levitator surface for cfMesh (exported as STL)
// Simple CSG: outer cylinder minus two shorter cylinders and a sphere.
SetFactory("OpenCASCADE");

// frequency
f = 25250;
// sound speed
c = 343;
// wavelength
l = c / f;
// distance
D = 0.00744; // 1.586*l;
// mesh size
N = 50;
// characteristic length (unused but kept from original)
gs = l / N;

// z-coordinate of sphere center (axis is y)
s = 0.004;
// sphere radius
r = 0.001;

// geometry
rT = 0.01;  // radius of transducer
rR = 0.019; // radius of reflector
bu = 0.01;  // air buffer
hS = 0.01;  // height of reflector/transducer stand

// Derived sizes
Rdom = rR + bu;
yMin = -hS;
yMax = D + hS;
Hdom = yMax - yMin;

// --- CSG solids ---
// Outer air domain (cylinder along +y)
volDom = newv; Cylinder(volDom) = {0, yMin, 0, 0, Hdom, 0, Rdom};

// Reflector cavity (bottom short cylinder)
volLev = newv; Cylinder(volLev) = {0, yMin, 0, 0, hS, 0, rR};

// Transducer cavity (top short cylinder)
volRef = newv; Cylinder(volRef) = {0, D, 0, 0, hS, 0, rT};

// Drop cavity (sphere on axis)
volDrop = newv; Sphere(volDrop) = {0, s, 0, r};

// Subtract cavities from the outer domain
volAir[] = BooleanDifference{ Volume{volDom}; Delete; }{ Volume{volLev, volRef, volDrop}; Delete; };

// Optional cleanup after booleans
Coherence;

transducer1[] = {9};
transducer2[] = {7};
openAir[] = {11,12,13};
reflector1[] = {5};
reflector2[] = {4};
dropWall[] = {10};

Physical Surface("transducer1") = {9};
Physical Surface("transducer2") = {7};
Physical Surface("openAir") = {11,12,13};
Physical Surface("reflector1") = {5};
Physical Surface("reflector2") = {4};
Physical Surface("dropWall") = {10};

Physical Volume("internal") = {1};

Mesh.CharacteristicLengthMin = gs;
Mesh.CharacteristicLengthMax = gs;

Mesh.ElementOrder = 1;
