// Gmsh project created on Fri Dec 20 10:38:20 2024
SetFactory("OpenCASCADE");

// Impact radius
ImR= 0.073;

// Droplet diameter
D = 0.002;

Box(1) = {ImR-5*D, -5*D, 0, 10*D, 10*D, 5*D};
//+
Physical Surface("Bottom", 13) = {5};
//+
Physical Surface("Mantle", 14) = {1, 4, 2, 3};
//+
Physical Surface("Top", 15) = {6};

Mesh.StlOneSolidPerSurface=2;

Save "squareSubstrate.stl";
