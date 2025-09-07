// vim: syntax=gmsh

SetFactory("Built-in");

cm = 1e-2;

Theta = Pi / 6;

lc = 1 * cm;

Point(1) = { 5 * cm, 0, -40 * cm};
Point(2) = {10 * cm, 0, -40 * cm};
Point(3) = {15 * cm, 0,   0 * cm};
Point(4) = { 5 * cm, 0,   0 * cm};

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

Extrude {-5 * cm, 0, 0}
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

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 2 - Theta}
  {
  Surface{44, 61, 83, 100, 117};
  Layers{(Pi / 2 - Theta) * 15 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Pi / 2}
  {
  Surface{139, 156, 178, 195, 212};
  Layers{(Pi / 2) * 15 * cm / lc};
  Recombine;
  }

Extrude {{0, 0, 1}, {0, 0, 0}, Pi}
  {
  Surface{234, 251, 273, 290, 307};
  Layers{Pi * 15 * cm / lc};
  Recombine;
  }

Physical Surface("inlet") = {177};

Physical Surface("outlet") = {112, 207, 302, 389};

Physical Surface("wall") = {31, 35, 56, 74, 78, 82, 116, 126, 130, 151, 169, 173, 211, 221, 225, 246,
                            264, 268, 272, 306, 316, 320, 340, 357, 361, 365, 393};

Physical Volume("domain") = {1:20};
