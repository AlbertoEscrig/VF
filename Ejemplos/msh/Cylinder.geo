
SetFactory("OpenCASCADE");

Point(1) = {0.5, 0, 0};

Extrude {0.5, 0, 0} { Point{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{1}; }

Point(5) = {2, 0, 0};
Point(6) = {2, Sin(Pi / 4), 0};

Curve(5) = {2, 5};
Curve(6) = {5, 6};
Curve(7) = {6, 4};

Curve Loop(7) = {5, 6, 7, -3};

Plane Surface(7) = {7};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{4}; }

Point(9)  = {0, 2, 0};
Point(10) = {Sin(Pi / 4), 2, 0};
Point(11) = {2, 2, 0};

Curve(11) = {4, 10};
Curve(12) = {10, 9};
Curve(13) = {9, 8};
Curve(14) = {6, 11};
Curve(15) = {11, 10};

Curve Loop(9) = {9, -13, -12, -11};
Plane Surface(9) = {9};

Curve Loop(11) = {7, 11, -15, -14};
Plane Surface(10) = {11};

Symmetry {1, 0, 0, 0} { Duplicata{ Surface{:}; } }

Coherence;

Transfinite Curve{:} = 20;
Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("cylinder") = {2, 16, 21, 28};
Physical Curve("down")     = {1, 5, 23, 25};
Physical Curve("up")       = {15, 19, 31, 32};
Physical Curve("right")    = {6, 14};
Physical Curve("left")     = {26, 33};

Physical Surface("domain") = Surface{:};
