
SetFactory("OpenCASCADE");

R = 1e-3;

Point(1) = {R, 0, 0};

Extrude {5 * R, 0, 0} { Point{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{1}; }

Point(5) = {20 * R, 0, 0};
Point(6) = {20 * R, 6 * R * Sin(Pi / 4), 0};

Curve(5) = {2, 5};
Curve(6) = {5, 6};
Curve(7) = {6, 4};

Curve Loop(7) = {5, 6, 7, -3};
Plane Surface(7) = {7};

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 4} { Curve{4}; }

Point(9)  = {0, 20 * R, 0};
Point(10) = {6 * R * Sin(Pi / 4), 20 * R, 0};
Point(11) = {20 * R, 20 * R, 0};

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
Extrude {80 * R, 0, 0} { Curve{6, 14}; }
Symmetry {0, 1, 0, 0} { Duplicata{ Surface{:}; } }

Coherence;

Transfinite Curve{7, 11, 14, 15, 60, 64, 70, 75, 77, 80, 82, 83, 90, 94, 96:98, 103, 106, 108, 109, 113} = 30;
Transfinite Curve{61:62, 65:67, 69, 71, 72, 76, 78, 79, 81, 85:87, 89, 91, 92, 95, 99, 100, 102, 104, 105, 107, 111} = 15;
Transfinite Curve{4, 63, 68, 73, 74, 88, 93, 101} = 30 Using Progression 1.05;
Transfinite Curve{57, 59, 84, 110, 112} = 150;

Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("cylinder") = {61, 66, 71, 78, 86, 91, 99, 104};
Physical Curve("inlet")    = {76, 83, 102, 109};
Physical Curve("outlet")   = {60, 85, 111, 113};
Physical Curve("up")       = {15, 59, 69, 81, 82};
Physical Curve("down")     = {95, 97, 107, 108, 112};

Physical Surface("domain") = Surface{:};
