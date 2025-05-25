// vim: syntax=gmsh

SetFactory("OpenCASCADE");

mm = 1e-3;

Point(1) = {0, 0, 0};

out[] = Extrude {100 * mm,       0, 0} { Point{1}; Layers{120}; };
out[] = Extrude {       0, 50 * mm, 0} { Curve{out[1]}; Layers{80};  Recombine; };
out[] = Extrude {       0, 10 * mm, 0} { Curve{out[0]}; Layers{20};  Recombine; };
        Extrude {50 * mm,        0, 0} { Curve{out[3]}; Layers{100}; Recombine; }
out[] = Extrude {      0, 240 * mm, 0} { Curve{out[0]}; Layers{380}; Recombine; };

Physical Curve("axis")   = {2, 5, 11};
Physical Curve("wall")   = {1, 3, 8, 9, 12, 13};
Physical Curve("nozzle") = {10};

Physical Surface("domain") = {1:4};
