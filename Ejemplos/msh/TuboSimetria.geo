// vim: syntax=gmsh

SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};

out[] = Extrude {0, 2.5e-2, 0} { Point{1}; Layers{25}; };
out[] = Extrude {10e-2, 0, 0} { Curve{out[1]}; Layers{100}; Recombine; };

S1 = out[1];

out[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi / 2 / 10} { Surface{S1}; Layers{1}; Recombine; };

S2 = out[0];
V1 = out[1];
S3 = out[2];
S4 = out[3];
S5 = out[4];

Physical Surface("symmetry") = {S1, S2};
Physical Surface("wall")     = {S3};
Physical Surface("inlet")    = {S4};
Physical Surface("outlet")   = {S5};

Physical Volume("domain") = {V1};
