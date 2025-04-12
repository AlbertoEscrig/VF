
SetFactory("OpenCASCADE");

p1 = newp; Point(p1) = {0, 0, 0};

out[] = Extrude {0, 2.5e-2, 0} { Point{p1}; Layers{25}; };

out[] = Extrude {10e-2, 0, 0} { Curve{out[1]}; Layers{100}; Recombine; };

xyplane = out[1];

out[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi / 2} { Surface{xyplane}; Layers{10}; Recombine; };

xzplane = out[0];
domain  = out[1];
wall    = out[2];
inlet   = out[3];
outlet  = out[4];

Mesh 3;

Periodic Surface{outlet}  = {inlet}   Translate {10e-2, 0, 0};
Periodic Surface{xzplane} = {xyplane} Rotate    {{1, 0, 0}, {0, 0, 0}, Pi / 2};

Physical Surface("wall")     = {wall};
Physical Surface("inlet")    = {inlet};
Physical Surface("outlet")   = {outlet};
Physical Surface("xy plane") = {xyplane};
Physical Surface("xz plane") = {xzplane};

Physical Volume("domain") = {domain};

