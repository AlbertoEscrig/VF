// vim: syntax=gmsh

SetFactory("OpenCASCADE");

N = 100;

Point(1) = {-1, -1, -1};

out[] = Extrude {2, 0, 0} { Point{1}; Layers{N}; };
out[] = Extrude {0, 2, 0} { Curve{out[1]}; Layers{N}; Recombine; };
out[] = Extrude {0, 0, 2} { Surface{out[1]}; Layers{N}; Recombine; };

Mesh 3;

S1 = Surface In BoundingBox {-1.01, -1.01, -1.01,  1.01,  1.01, -0.99};
S2 = Surface In BoundingBox {-1.01, -1.01,  0.99,  1.01,  1.01,  1.01};
S3 = Surface In BoundingBox {-1.01, -1.01, -1.01,  1.01, -0.99,  1.01};
S4 = Surface In BoundingBox {-1.01,  0.99, -1.01,  1.01,  1.01,  1.01};
S5 = Surface In BoundingBox {-1.01, -1.01, -1.01, -0.99,  1.01,  1.01};
S6 = Surface In BoundingBox { 0.99, -1.01, -1.01,  1.01,  1.01,  1.01};

Periodic Surface{S2} = {S1} Translate {0, 0, 2};
Periodic Surface{S4} = {S3} Translate {0, 2, 0};
Periodic Surface{S6} = {S5} Translate {2, 0, 0};

Physical Surface("xy plane") = {S1, S2};
Physical Surface("xz plane") = {S3, S4};
Physical Surface("yz plane") = {S5, S6};

Physical Volume("domain") = Volume{:};
