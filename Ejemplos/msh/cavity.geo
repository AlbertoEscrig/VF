
SetFactory("OpenCASCADE");

p1 = newp; Point(p1) = {0, 0, 0};

out[] = Extrude {0.1, 0, 0} { Point{p1}; Layers{20}; };

out[] = Extrude {0, 0.1, 0} { Curve{out[1]}; Layers{20}; Recombine; };

Physical Curve("fixedWalls") = {1, 2, 3};
Physical Curve("movingWall") = {4};

Physical Surface("domain") = {1};
