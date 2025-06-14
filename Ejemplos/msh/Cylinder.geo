// vim: syntax=gmsh

SetFactory("Built-in");

Point(1) = {0.5, 0, 0};

Extrude {0.5, 0, 0} { Point{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{1}; }

Point(5) = {2, 0, 0};
Point(6) = {2, Sin(Pi / 4), 0};

Curve(5) = {5, 2};
Curve(6) = {6, 5};
Curve(7) = {4, 6};

Curve Loop(7) = {4:7};
Plane Surface(7) = {7};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{2}; }

Point(12) = {0, 2, 0};
Point(13) = {Sin(Pi / 4), 2, 0};
Point(14) = {2, 2, 0};

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

Coherence;

Transfinite Curve{:} = 20;
Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("cylinder") = {3, 9, 20, 30};
Physical Curve("down")     = {1, 5, 17, 25};
Physical Curve("up")       = {12, 15, 34, 39};
Physical Curve("right")    = {6, 14};
Physical Curve("left")     = {24, 38};

Physical Surface("domain") = Surface{:};
