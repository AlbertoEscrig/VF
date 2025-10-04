// vim: syntax=gmsh

SetFactory("OpenCASCADE");

mm = 1e-3;

Point(1) = {0, 0, 0};

C1[] = {};
C2[] = {};

out[] = Extrude {100 * mm, 0, 0} { Point{1}; Layers{120}; };
C1[] += out[1];

out[] = Extrude {0, 50 * mm, 0} { Curve{out[1]}; Layers{80};  Recombine; };
C1[] += out[3];
C2[] += out[2];

out[] = Extrude {0, 10 * mm, 0} { Curve{out[0]}; Layers{20};  Recombine; };
C2[] += out[2];
C3 = out[3];

out[] = Extrude {0, 240 * mm, 0} { Curve{out[0]}; Layers{380}; Recombine; };
C1[] += {out[0], out[3]};
C2[] += out[2];

out[] = Extrude {50 * mm, 0, 0} { Curve{C3}; Layers{100}; Recombine; };
C1[] += {out[2], out[3]};

Physical Curve("wall")   = C1[];
Physical Curve("axis")   = C2[];
Physical Curve("nozzle") = out[0];

Physical Surface("domain") = Surface{:};
