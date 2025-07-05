// vim: syntax=gmsh

SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};

C[] = {};

out[] = Extrude {0.45, 0, 0} { Point{1};      Layers{45}; }; C[] += out[1];
out[] = Extrude {0.10, 0, 0} { Point{out[0]}; Layers{10}; }; C[] += out[1];
out[] = Extrude {0.45, 0, 0} { Point{out[0]}; Layers{45}; }; C[] += out[1];

out[] = Extrude {0, 0.45, 0} { Curve{C[]}; Layers{45}; Recombine; };

C[] += {out[2], out[11]};
S1[] = {out[1], out[5], out[9]};

out[] = Extrude {0, 0.10, 0} { Curve{out[0], out[4], out[8]}; Layers{10}; Recombine; };

C[] += {out[2], out[11]};
S1[] += {out[1], out[9]};
S2 = out[5];

out[] = Extrude {0, 0.45, 0} { Curve{out[0], out[4], out[8]}; Layers{45}; Recombine; };

C[] += {out[0], out[2], out[4], out[8], out[11]};
S1[] += {out[1], out[5], out[9]};

Physical Curve("boundary") = C[];
Physical Surface("domain") = S1[];
Physical Surface("seed")   = S2;
