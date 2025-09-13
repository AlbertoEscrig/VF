// vim: syntax=gmsh

SetFactory("Built-in");

cm = 1e-2;

Theta = Pi / 3;

lc = 1 * cm;

Point(1) = {0,  5 * cm, -40 * cm};
Point(2) = {0, 10 * cm, -40 * cm};
Point(3) = {0, 15 * cm,   0 * cm};
Point(4) = {0,  5 * cm,   0 * cm};

Curve(1) = {1, 2};
Curve(2) = {2, 3};
Curve(3) = {3, 4};
Curve(4) = {4, 1};

Transfinite Curve{1, 3} = 10 * cm / lc;
Transfinite Curve{2, 4} = 40 * cm / lc;

Curve Loop(1) = {1:4};
Plane Surface(1) = {1};

Transfinite Surface{1};
Recombine Surface{1};

Extrude {0, -5 * cm, 0}
  {
  Curve{4};
  Layers{5 * cm / lc};
  Recombine;
  }

Extrude {0, 0, 10 * cm}
  {
  Curve{3, 6};
  Layers{10 * cm / lc};
  Recombine;
  }

Extrude {0, 0, 10 * cm}
  {
  Curve{13};
  Layers{10 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Theta}
  {
  Surface{1, 8, 12, 16, 20};
  Layers{Theta * 15 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Pi - Theta / 2}
  {
  Surface{42, 59, 81, 98, 115};
  Layers{(Pi - Theta / 2) * 15 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Pi - Theta / 2}
  {
  Surface{137, 154, 176, 193, 210};
  Layers{(Pi - Theta / 2) * 15 * cm / lc};
  Recombine;
  }

Physical Surface("inlet")  = {80};
Physical Surface("outlet") = {110, 205, 292};
Physical Surface("wall")   = {29, 33, 54, 72, 76, 114, 124, 128, 149, 167, 171,
                              175, 209, 219, 223, 243, 260, 264, 268, 296};

Physical Volume("domain") = Volume{:};
