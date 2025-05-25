// vim: syntax=gmsh

SetFactory("OpenCASCADE");

p1 = newp; Point(p1) = {0, 0, 0};

out[] = Extrude {0.1, 0, 0} { Point{p1}; Layers{20}; };

L1 = out[1];

out[] = Extrude {0, 0.1, 0} { Curve{L1}; Layers{20}; Recombine; };

L2 = out[0];
S1 = out[1];
L3 = out[2];
L4 = out[3];

Physical Curve("fixedWalls") = {L1, L3, L4};
Physical Curve("movingWall") = {L2};

Physical Surface("domain") = {S1};
