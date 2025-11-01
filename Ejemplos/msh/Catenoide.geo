// vim: syntax=gmsh

SetFactory("OpenCASCADE");

Point(1) = {  0,   0, 0};
Point(2) = {0.5,   0, 0};
Point(3) = {  1,   0, 0};
Point(4) = {  0, 0.5, 0};
Point(5) = {  0,   1, 0};

Curve(1)  = {2, 3};
Circle(2) = {3, 1, 5};
Curve(3)  = {5, 4};
Circle(4) = {4, 1, 2};

Curve Loop(1) = {1:4};
Plane Surface(1) = {1};

Transfinite Curve{1, 3} = 50;
Transfinite Curve{2, 4} = 200;

Transfinite Surface{1};
Recombine Surface{1};

Physical Curve("symmetry") = {1, 3};
Physical Curve("outer")    = {2};
Physical Curve("inner")    = {4};

Physical Surface("domain") = {1};
