
SetFactory("OpenCASCADE");

p1 = newp; Point(p1) = {0, 0, 0};

out[] = Extrude {0.1, 0, 0} { Point{p1}; };

out[] = Extrude {0, 0.1, 0} { Curve{out[1]}; };

Transfinite Curve{:} = 75 Using Bump 0.2;

Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("adiabaticWall") = {1, 4};
Physical Curve("hotWall") = {2};
Physical Curve("coldWall") = {3};

Physical Surface("domain") = {1};
