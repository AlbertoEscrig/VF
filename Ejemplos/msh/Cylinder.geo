// vim: syntax=gmsh

SetFactory("OpenCASCADE");

R = 0.5;
L = 2;

Point(1)  = { 0, 0, 0};
Point(2)  = {-L, 0, 0};
Point(3)  = {-R, 0, 0};
Point(4)  = { 0, R, 0};
Point(5)  = { R, 0, 0};
Point(6)  = { L, 0, 0};
Point(7)  = { L, L, 0};
Point(8)  = { R, L, 0};
Point(9)  = { 0, L, 0};
Point(10) = {-R, L, 0};
Point(11) = {-L, L, 0};

Curve(1)  = {2, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Curve(4)  = {5, 6};
Curve(5)  = {6, 7};
Curve(6)  = {7, 8};
Curve(7)  = {8, 9};
Curve(8)  = {9, 10};
Curve(9)  = {10, 11};
Curve(10) = {11, 2};

Curve Loop(1) = {1:10};
Plane Surface(1) = {1};

Transfinite Curve{1, 4, 5, 6, 9, 10} = 75;
Transfinite Curve{2, 3, 7, 8} = 25;

Transfinite Surface{1} = {2, 6, 7, 11};
Recombine Surface{1};

Mesh.Smoothing = 100;

Physical Curve("cylinder") = {2, 3};
Physical Curve("down")     = {1, 4};
Physical Curve("up")       = {6, 7, 8, 9};
Physical Curve("right")    = {5};
Physical Curve("left")     = {10};

Physical Surface("domain") = {1};
