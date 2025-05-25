// vim: syntax=gmsh

SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};

out[] = Extrude {0.6,   0, 0} { Point{1}; Layers{90}; };
out[] = Extrude {  0, 0.2, 0} { Curve{out[1]}; Layers{30}; Recombine; };
out[] = Extrude {  0, 0.8, 0} { Curve{out[0]}; Layers{120}; Recombine; };
out[] = Extrude {2.4,   0, 0} { Curve{out[3]}; Layers{360}; Recombine; };

Physical Curve("inlet")    = {2, 5};
Physical Curve("outlet")   = {10};
Physical Curve("bottom")   = {1};
Physical Curve("obstacle") = {3, 8};
Physical Curve("top")      = {7, 9};

Physical Surface("domain") = {1:3};
