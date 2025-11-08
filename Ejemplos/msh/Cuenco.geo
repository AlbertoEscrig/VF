// vim: syntax=gmsh

SetFactory("OpenCASCADE");

R = 0.1;
alfa = Pi / 4;

c = R * Cos(alfa);
s = R * Sin(alfa);

Point(1) = { 0, 0, 0};
Point(2) = { R, 0, 0};
Point(3) = { c, s, 0};
Point(4) = {-c, s, 0};
Point(5) = {-R, 0, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Curve(4)  = {5, 2};

Curve Loop(1) = {1:4};
Plane Surface(1) = {1};

Mesh.Algorithm = 11;
Mesh.MeshSizeMax = 2e-3;

Physical Curve("cold")      = 1;
Physical Curve("adiabatic") = 2;
Physical Curve("hot")       = 3;
Physical Curve("symmetry")  = 4;

Physical Surface("domain") = 1;
