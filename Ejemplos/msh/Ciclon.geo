
cm = 1e-2;

Theta = Pi / 6;

lc = 1 * cm;

p1 = newp; Point(p1) = { 5 * cm, 0, -40 * cm};
p2 = newp; Point(p2) = {10 * cm, 0, -40 * cm};
p3 = newp; Point(p3) = {15 * cm, 0,   0 * cm};
p4 = newp; Point(p4) = { 5 * cm, 0,   0 * cm};

c1 = newc; Curve(c1) = {p1, p2};
c2 = newc; Curve(c2) = {p2, p3};
c3 = newc; Curve(c3) = {p3, p4};
c4 = newc; Curve(c4) = {p4, p1};

Transfinite Curve{c1, c3} = 10 * cm / lc;
Transfinite Curve{c2, c4} = 40 * cm / lc;

cc = newcl; Curve Loop(cc) = {c1, c2, c3, c4};

s1 = news; Plane Surface(s1) = {cc};

Transfinite Surface{s1};
Recombine Surface{s1};

Extrude {-5 * cm, 0, 0}
  {
  Curve{c4};
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

