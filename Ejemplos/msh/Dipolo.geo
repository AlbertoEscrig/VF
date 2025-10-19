// vim: syntax=gmsh

SetFactory("OpenCASCADE");

R = 0.2;
L = 1.0;
H = 5.0;

lcF = 0.05;
lcC = 0.2;

Circle(1) = { L / 2, 0, 0, R};
Circle(2) = {-L / 2, 0, 0, R};
Circle(3) = {     0, 0, 0, H};

Curve Loop(1) = {1};
Curve Loop(2) = {2};
Curve Loop(3) = {3};

Plane Surface(1) = {1, 2, 3};

Field[1] = Distance;
Field[1].CurvesList = {1, 2};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lcF;
Field[2].LcMax = lcC;
Field[2].DistMin = R;
Field[2].DistMax = H;

Background Field = 2;

Physical Curve("positive")  = {1};
Physical Curve("negative")  = {2};
Physical Curve("far field") = {3};

Physical Surface("domain") = {1};
