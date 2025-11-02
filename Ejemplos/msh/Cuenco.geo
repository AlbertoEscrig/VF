// vim: syntax=gmsh

SetFactory("Built-in");

Point(1) = {                 0,                 0, 0};
Point(2) = {               0.1,                 0, 0};
Point(3) = { 0.1 * Cos(Pi / 4), 0.1 * Sin(Pi / 4), 0};
Point(4) = {-0.1 * Cos(Pi / 4), 0.1 * Sin(Pi / 4), 0};
Point(5) = {              -0.1,                 0, 0};

Curve(1)  = {1, 2};
Circle(2) = {2, 1, 3};
Curve(3)  = {3, 1};
Circle(4) = {3, 1, 4};
Curve(5)  = {4, 1};
Circle(6) = {4, 1, 5};
Curve(7)  = {5, 1};

Curve Loop(1) = {1, 2, 3};  Plane Surface(1) = {1};
Curve Loop(2) = {-3, 4, 5}; Plane Surface(2) = {2};
Curve Loop(3) = {-5, 6, 7}; Plane Surface(3) = {3};

Transfinite Curve{1:3, 5:7} = 100;
Transfinite Curve{4} = 200;

Transfinite Surface{1} = {1, 2, 3};
Transfinite Surface{2} = {1, 3, 4};
Transfinite Surface{3} = {1, 4, 5};

Recombine Surface{1:3};

Physical Curve("symmetry")  = {1, 7};
Physical Curve("cold")      = {2};
Physical Curve("adiabatic") = {4};
Physical Curve("hot")       = {6};

Physical Surface("domain") = {1:3};
