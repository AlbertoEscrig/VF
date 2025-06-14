// vim: syntax=gmsh

SetFactory("Built-in");

R = 1e-3;

Point(1) = {R, 0, 0};

Extrude {5 * R, 0, 0} { Point{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{1}; }

Point(5) = {20 * R, 0, 0};
Point(6) = {20 * R, 6 * R * Sin(Pi / 4), 0};

Curve(5) = {5, 2};
Curve(6) = {6, 5};
Curve(7) = {4, 6};

Curve Loop(7) = {4:7};
Plane Surface(7) = {7};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{2}; }

Point(12)  = {0, 20 * R, 0};
Point(13) = {6 * R * Sin(Pi / 4), 20 * R, 0};
Point(14) = {20 * R, 20 * R, 0};

Curve(11) = {13, 4};
Curve(12) = {12, 13};
Curve(13) = {10, 12};
Curve(14) = {6, 14};
Curve(15) = {14, 13};

Curve Loop(12) = {10:13};
Plane Surface(12) = {12};

Curve Loop(13) = {7, 11, 14, 15};
Plane Surface(13) = {13};

Symmetry {1, 0, 0, 0} { Duplicata{ Surface{:}; } }
Extrude {30 * R, 0, 0} { Curve{6, 14}; }
Symmetry {0, 1, 0, 0} { Duplicata{ Surface{:}; } }

Coherence;

Transfinite Curve{11, 13, 14, 23, 25, 35, 38, 39, 44, 65, 67, 70, 80, 92, 95, 96, 106} = 30;
Transfinite Curve{3, 4, 6, 9, 10, 12, 18, 20, 24, 28, 30, 34, 40, 50, 52, 56, 60, 62, 66, 75, 77, 81, 85, 87, 91, 101} = 15;
Transfinite Curve{1, 2, 8, 17, -19, -51, -61, -76} = 30 Using Progression 1.05;
Transfinite Curve{-5, 7, -15, 55, -71} = 30 Using Progression 1.02;
Transfinite Curve{41, 42, 46, 102, 105} = 40;

Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("cylinder") = {3, 9, 20, 30, 52, 62, 77, 87};
Physical Curve("inlet")    = {24, 38, 81, 95};
Physical Curve("outlet")   = {40, 44, 101, 106};
Physical Curve("up")       = {12, 15, 34, 39, 46};
Physical Curve("down")     = {66, 71, 91, 96, 105};

Physical Surface("domain") = Surface{:};
