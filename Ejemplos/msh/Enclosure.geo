// vim: syntax=gmsh

SetFactory("OpenCASCADE");

L = 1.0;
N = 60;

Point(1) = {0, 0, 0};

out[] = Extrude {L / 3, 0, 0} { Point{1}; Layers{N / 3}; };
C[] = out[1];

out[] = Extrude {L / 3, 0, 0} { Point{out[0]}; Layers{N / 3}; };
C[] += out[1];

out[] = Extrude {L / 3, 0, 0} { Point{out[0]}; Layers{N / 3}; };
C[] += out[1];

out[] = Extrude {0, L / 3, 0} { Curve{C[]}; Layers{N / 3}; Recombine; };
S[] = {out[1], out[5], out[9]};

out[] = Extrude {0, L / 3, 0} { Curve{out[0], out[4], out[8]}; Layers{N /3}; Recombine; };
S[] += {out[1], out[9]};

Delete { Surface{out[5]}; }

out[] = Extrude {0, L / 3, 0} { Curve{out[0], out[4], out[8]}; Layers{N / 3}; Recombine; };
S[] += {out[1], out[5], out[9]};

out[] = Extrude {0, 0, L} { Surface{S[]}; Layers{N}; Recombine; };

e = 1e-4;

Physical Surface("hot") = {Surface In BoundingBox {-e, -e, -e, e, L + e, L + e}};

Physical Surface("cold") =
{
  Surface In BoundingBox {          - e,           - e,   - e,     L     + e,     L     + e,     e},
  Surface In BoundingBox {          - e,           - e,   - e,     L     + e,             e, L + e},
  Surface In BoundingBox {          - e,           - e, L - e,     L     + e,     L     + e, L + e},
  Surface In BoundingBox {          - e,     L     - e,   - e,     L     + e,     L     + e, L + e},
  Surface In BoundingBox {    L / 3 - e,     L / 3 - e,   - e,     L / 3 + e, 2 * L / 3 + e, L + e},
  Surface In BoundingBox {    L / 3 - e,     L / 3 - e,   - e, 2 * L / 3 + e,     L / 3 + e, L + e},
  Surface In BoundingBox {    L / 3 - e, 2 * L / 3 - e,   - e, 2 * L / 3 + e, 2 * L / 3 + e, L + e},
  Surface In BoundingBox {2 * L / 3 - e,     L / 3 - e,   - e, 2 * L / 3 + e, 2 * L / 3 + e, L + e},
  Surface In BoundingBox {    L     - e,           - e,   - e,     L     + e,     L     + e, L + e}
};

Physical Volume("domain") = Volume{:};
