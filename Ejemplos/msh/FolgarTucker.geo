// vim: syntax=gmsh

SetFactory("OpenCASCADE");

H = 0.05;
L = 0.4;

Point(1) = {0, -H, 0};
Point(2) = {L, -H, 0};
Point(3) = {L,  H, 0};
Point(4) = {0,  H, 0};

Curve(1) = {1, 2};
Curve(2) = {2, 3};
Curve(3) = {3, 4};
Curve(4) = {4, 1};

Curve Loop(1) = {1:4};
Plane Surface(1) = {1};

Transfinite Curve{1, 3} = 200;
Transfinite Curve{2, 4} = 100 Using Bump 0.5;

Transfinite Surface{1};
Recombine Surface{1};

Physical Curve("wall")   = {1, 3};
Physical Curve("outlet") = 2;
Physical Curve("inlet")  = 4;

Physical Surface("domain") = 1;
