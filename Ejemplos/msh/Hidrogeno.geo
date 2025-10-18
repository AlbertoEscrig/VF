// vim: syntax=gmsh

SetFactory("OpenCASCADE");

L  = 30;
Lz = 1.1 * L;
N  = 70;

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, L, 0};
Point(4) = {0, L, 0};
Point(5) = {0, 0, Lz};
Point(6) = {L, 0, Lz};
Point(7) = {L, L, Lz};
Point(8) = {0, L, Lz};

Curve(1) = {1, 2}; Curve(2)  = {1, 4}; Curve(3)  = {2, 3}; Curve(4)  = {4, 3};
Curve(5) = {1, 5}; Curve(6)  = {2, 6}; Curve(7)  = {4, 8}; Curve(8)  = {3, 7};
Curve(9) = {5, 6}; Curve(10) = {5, 8}; Curve(11) = {6, 7}; Curve(12) = {8, 7};

Curve Loop(1) = {2,  7, -10,  -5}; Plane Surface(1) = {1};
Curve Loop(2) = {3,  8, -11,  -6}; Plane Surface(2) = {2};
Curve Loop(3) = {1,  6,  -9,  -5}; Plane Surface(3) = {3};
Curve Loop(4) = {4,  8, -12,  -7}; Plane Surface(4) = {4};
Curve Loop(5) = {1,  3,  -4,  -2}; Plane Surface(5) = {5};
Curve Loop(6) = {9, 11, -12, -10}; Plane Surface(6) = {6};

Surface Loop(1) = {1:6}; Volume(1) = {1};

Transfinite Curve{1:12} = N Using Progression 1.01;

Transfinite Surface{1:6};
Recombine Surface{1:6};

Transfinite Volume{1};

Physical Surface("symmetry")  = {1, 3, 5};
Physical Surface("far field") = {2, 4, 6};

Physical Volume("domain") = {1};
