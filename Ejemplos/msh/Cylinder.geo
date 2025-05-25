// vim: syntax=gmsh

SetFactory("Built-in");

Point(1) = {0.5, 0, 0};

Extrude {0.5, 0, 0} { Point{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{1}; }

Point(5) = {2, 0, 0};
Point(6) = {2, Sin(Pi / 4), 0};

Curve(5) = {2, 5};
Curve(6) = {5, 6};
Curve(7) = {6, 4};

Curve Loop(7) = {-4, 5, 6, 7};
Plane Surface(7) = {7};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{2}; }

Point(12) = {0, 2, 0};
Point(13) = {Sin(Pi / 4), 2, 0};
Point(14) = {2, 2, 0};

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

Coherence;

Transfinite Curve{:} = 20;
Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("cylinder") = {3, 9, 20, 30};
Physical Curve("down")     = {1, 5, 17, 23};
Physical Curve("up")       = {12, 15, 34, 39};
Physical Curve("right")    = {6, 14};
Physical Curve("left")     = {24, 40};

Physical Surface("domain") = Surface{:};
