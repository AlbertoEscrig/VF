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

Curve Loop(5) = {1:4};

Plane Surface(6) = {5};

Transfinite Surface{6};
Recombine Surface{6};

Extrude {0, -5 * cm, 0}
  {
  Curve{4};
  Layers{5 * cm / lc};
  Recombine;
  }

Extrude {0, 0, 10 * cm}
  {
  Curve{3, 8};
  Layers{10 * cm / lc};
  Recombine;
  }

Extrude {0, 0, 10 * cm}
  {
  Curve{15};
  Layers{10 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Theta}
  {
  Surface{6, 10, 14, 18, 22};
  Layers{Theta * 15 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Pi - Theta / 2}
  {
  Surface{44, 61, 83, 100, 117};
  Layers{(Pi - Theta / 2) * 15 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Pi - Theta / 2}
  {
  Surface{139, 156, 178, 195, 212};
  Layers{(Pi - Theta / 2) * 15 * cm / lc};
  Recombine;
  }

Physical Surface("inlet")  = {82};
Physical Surface("outlet") = {112, 207, 294};
Physical Surface("wall")   = {31, 35, 56, 74, 78, 116, 126, 130, 151, 169, 173,
                              177, 211, 221, 225, 245, 262, 266, 270, 298};

Physical Volume("domain") = Volume{:};
