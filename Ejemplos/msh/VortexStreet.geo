// vim: syntax=gmsh

SetFactory("OpenCASCADE");

R = 1e-3;

Point(1) = {0, 0, 0};

Point(2) = { R * Sqrt(2) / 2,  R * Sqrt(2) / 2, 0};
Point(3) = {-R * Sqrt(2) / 2,  R * Sqrt(2) / 2, 0};
Point(4) = {-R * Sqrt(2) / 2, -R * Sqrt(2) / 2, 0};
Point(5) = { R * Sqrt(2) / 2, -R * Sqrt(2) / 2, 0};

Point(6) = { 3 * R * Sqrt(2),  3 * R * Sqrt(2), 0};
Point(7) = {-3 * R * Sqrt(2),  3 * R * Sqrt(2), 0};
Point(8) = {-3 * R * Sqrt(2), -3 * R * Sqrt(2), 0};
Point(9) = { 3 * R * Sqrt(2), -3 * R * Sqrt(2), 0};

Point(10) = {-20 * R, 3 * R * Sqrt(2), 0};
Point(11) = { 20 * R, 3 * R * Sqrt(2), 0};
Point(12) = { 50 * R, 3 * R * Sqrt(2), 0};

Point(13) = {-20 * R, -3 * R * Sqrt(2), 0};
Point(14) = { 20 * R, -3 * R * Sqrt(2), 0};
Point(15) = { 50 * R, -3 * R * Sqrt(2), 0};

Point(16) = {-20 * R,           20 * R, 0};
Point(17) = { -3 * R * Sqrt(2), 20 * R, 0};
Point(18) = {  3 * R * Sqrt(2), 20 * R, 0};
Point(19) = { 20 * R,           20 * R, 0};
Point(20) = { 50 * R,           20 * R, 0};

Point(21) = {-20 * R,           -20 * R, 0};
Point(22) = { -3 * R * Sqrt(2), -20 * R, 0};
Point(23) = {  3 * R * Sqrt(2), -20 * R, 0};
Point(24) = { 20 * R,           -20 * R, 0};
Point(25) = { 50 * R,           -20 * R, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Curve(9)  = {2, 6};
Curve(10) = {3, 7};
Curve(11) = {4, 8};
Curve(12) = {5, 9};

Curve(13) = {10, 7};
Curve(14) = {6, 11};
Curve(15) = {11, 12};

Curve(16) = {13, 8};
Curve(17) = {9, 14};
Curve(18) = {14, 15};

Curve(19) = {16, 17};
Curve(20) = {17, 18};
Curve(21) = {18, 19};
Curve(22) = {19, 20};

Curve(23) = {21, 22};
Curve(24) = {22, 23};
Curve(25) = {23, 24};
Curve(26) = {24, 25};

Curve(27) = {11, 14};

Curve(28) = {10, 13};
Curve(29) = {10, 16};
Curve(30) = {13, 21};

Curve(31) = {12, 15};
Curve(32) = {12, 20};
Curve(33) = {15, 25};

Curve Loop(1) = {1, 10, -5, -9};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 11, -6, -10};
Plane Surface(2) = {2};

Curve Loop(3) = {3, 12, -7, -11};
Plane Surface(3) = {3};

Curve Loop(4) = {4, 9, -8, -12};
Plane Surface(4) = {4};

Curve Loop(5) = {6, -16, -28, 13};
Plane Surface(5) = {5};

Curve Loop(6) = {8, 14, 27, -17};
Plane Surface(6) = {6};

Curve Loop(7) = {15, 31, -18, -27};
Plane Surface(7) = {7};

Curve Loop(8) = {5, -13, 29, 19, 20, 21, 22, -32, -15, -14};
Plane Surface(8) = {8};

Curve Loop(9) = {7, 17, 18, 33, -26, -25, -24, -23, -30, 16};
Plane Surface(9) = {9};

Transfinite Curve{1:8, 13, 16, 19, 20, 23, 24, 27:33} = 30;
Transfinite Curve{9:12} = 30 Using Progression 1.05;
Transfinite Curve{14, 17, 21, 25} = 30 Using Progression 1.02;
Transfinite Curve{15, 18, 22, 26} = 40;

Transfinite Surface{1:7};
Transfinite Surface{8} = {10, 12, 16, 20};
Transfinite Surface{9} = {13, 15, 21, 25};

Recombine Surface{:};

Physical Curve("cylinder") = {1:4};
Physical Curve("inlet")    = {28:30};
Physical Curve("outlet")   = {31:33};
Physical Curve("up")       = {19:22};
Physical Curve("down")     = {23:26};

Physical Surface("domain") = {1:9};
