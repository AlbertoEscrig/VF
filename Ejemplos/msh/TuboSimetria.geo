
cm = 1e-2;

Point(1) = {0, 0, 0};

out[] = Extrude {0, 2.5 * cm, 0} { Point{1}; Layers{25}; };

out[] = Extrude {10 * cm, 0, 0} { Curve{out[1]}; Layers{100}; Recombine; };

out[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi / 2 / 10} { Surface{out[1]}; Layers{1}; Recombine; };

Physical Surface("wall")     = {17};
Physical Surface("inlet")    = {13};
Physical Surface("outlet")   = {20};
Physical Surface("symmetry") = {5, 22};

Physical Volume("domain") = {1};
