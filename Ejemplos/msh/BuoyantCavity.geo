// vim: syntax=gmsh

SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};

out[] = Extrude {0.1, 0, 0} { Point{1}; };

L1 = out[1];

out[] = Extrude {0, 0.1, 0} { Curve{L1}; };

L2 = out[0];
S1 = out[1];
L3 = out[2];
L4 = out[3];

Transfinite Curve{:} = 75 Using Bump 0.2;

Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("adiabaticWall") = {L1, L2};
Physical Curve("hotWall")       = {L3};
Physical Curve("coldWall")      = {L4};

Physical Surface("domain") = {S1};
