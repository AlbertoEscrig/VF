// vim: syntax=gmsh

SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};

out[] = Extrude {0.6, 0, 0} { Point{1}; Layers{90}; };
C1 = out[1];

out[] = Extrude {0, 0.2, 0} { Curve{C1}; Layers{30}; Recombine; };
C2[] = out[2];
C3[] = out[3];

out[] = Extrude {0, 0.8, 0} { Curve{out[0]}; Layers{120}; Recombine; };
C2[] += out[2];
C4[] = out[0];

out[] = Extrude {2.4, 0, 0} { Curve{out[3]}; Layers{360}; Recombine; };
C3[] += out[2];
C4[] += out[3];
C5 = out[0];

Physical Curve("bottom")   = C1;
Physical Curve("inlet")    = C2[];
Physical Curve("obstacle") = C3[];
Physical Curve("top")      = C4[];
Physical Curve("outlet")   = C5;

Physical Surface("domain") = Surface{:};
