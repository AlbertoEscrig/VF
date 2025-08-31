// vim: syntax=gmsh

SetFactory("OpenCASCADE");

N = 100;

Point(1) = {-1, -1, -1};

out[] = Extrude {0, 2, 0} { Point{1}; Layers{N}; };

out[] = Extrude {2, 0, 0} { Curve{out[1]}; Layers{N}; Recombine; };

out[] = Extrude {0, 0, 2} { Surface{out[1]}; Layers{N}; Recombine; };

Mesh 3;

Periodic Surface{6} = {1} Translate {0, 0, 2};
Periodic Surface{3} = {2} Translate {0, 2, 0};
Periodic Surface{5} = {4} Translate {2, 0, 0};

Physical Surface("xy plane") = {1, 6};
Physical Surface("xz plane") = {2, 3};
Physical Surface("yz plane") = {4, 5};

Physical Volume("domain") = {1};
