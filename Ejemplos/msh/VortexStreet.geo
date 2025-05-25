// vim: syntax=gmsh

SetFactory("Built-in");

R = 1e-3;

Point(1) = {R, 0, 0};

Extrude {5 * R, 0, 0} { Point{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{1}; }

Point(5) = {20 * R, 0, 0};
Point(6) = {20 * R, 6 * R * Sin(Pi / 4), 0};

Curve(5) = {2, 5};
Curve(6) = {5, 6};
Curve(7) = {6, 4};

Curve Loop(7) = {-4, 5, 6, 7};
Plane Surface(7) = {7};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{2}; }

Point(12)  = {0, 20 * R, 0};
Point(13) = {6 * R * Sin(Pi / 4), 20 * R, 0};
Point(14) = {20 * R, 20 * R, 0};

Curve(11) = {4, 13};
Curve(12) = {12, 13};
Curve(13) = {10, 12};
Curve(14) = {14, 6};
Curve(15) = {13, 14};

Curve Loop(12) = {10, -11, 12, 13};
Plane Surface(12) = {12};

Curve Loop(13) = {7, 11, 14, 15};
Plane Surface(13) = {13};

Symmetry {1, 0, 0, 0} { Duplicata{ Surface{:}; } }
Extrude {80 * R, 0, 0} { Curve{6, 14}; }
Symmetry {0, 1, 0, 0} { Duplicata{ Surface{:}; } }

Coherence;

Transfinite Curve{5, 7, 11, 13:15, 23, 25, 35, 39, 40, 45, 58, 66, 68, 72, 73, 83, 93, 97, 98, 107} = 30;
Transfinite Curve{3, 4, 6, 9, 10, 12, 18, 20, 24, 28, 30, 34, 41, 51, 53, 57, 61, 63, 67, 76, 78, 82, 86, 88, 92, 102} = 15;
Transfinite Curve{1, 2, 8, 17, -19, -52, -62, -77} = 30 Using Progression 1.05;
Transfinite Curve{42, 43, 46, 101, 108} = 150;

Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("cylinder") = {3, 9, 20, 30, 53, 63, 78, 88};
Physical Curve("inlet")    = {24, 40, 82, 98};
Physical Curve("outlet")   = {41, 45, 102, 107};
Physical Curve("up")       = {12, 15, 34, 39, 46};
Physical Curve("down")     = {67, 72, 92, 97, 108};

Physical Surface("domain") = Surface{:};
