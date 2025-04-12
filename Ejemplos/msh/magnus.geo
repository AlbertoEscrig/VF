
SetFactory("OpenCASCADE");

R = 10e-2;

p1 = newp; Point(p1) = { R, 0, 0 };

out[] = Extrude {5, 0, 0} { Point{p1}; };

out[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi}  { Curve{out[1]}; };

out[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi}  { Curve{out[0]}; };

Coherence;

Transfinite Curve{4, 7} = 100 Using Progression 1.05;

Transfinite Curve{5, 6, 8, 9} = 100;

Transfinite Surface "*";

Recombine Surface "*";

Physical Curve("ball") = {8, 5};

Physical Curve("inlet") = {9};

Physical Curve("outlet") = {6};

Physical Surface("domain") = {1, 2};
