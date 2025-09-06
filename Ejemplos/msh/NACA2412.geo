// vim: syntax=gmsh

SetFactory("OpenCASCADE");

lcF = 0.005;
lcC = 0.2;

Point(1)  = {1.0000,  0.0013, 0, lcF};
Point(2)  = {0.9500,  0.0113, 0, lcF};
Point(3)  = {0.9000,  0.0209, 0, lcF};
Point(4)  = {0.8000,  0.0385, 0, lcF};
Point(5)  = {0.7000,  0.0543, 0, lcF};
Point(6)  = {0.6000,  0.0679, 0, lcF};
Point(7)  = {0.5000,  0.0788, 0, lcF};
Point(8)  = {0.4000,  0.0864, 0, lcF};
Point(9)  = {0.3000,  0.0891, 0, lcF};
Point(10) = {0.2500,  0.0876, 0, lcF};
Point(11) = {0.2000,  0.0840, 0, lcF};
Point(12) = {0.1500,  0.0777, 0, lcF};
Point(13) = {0.1000,  0.0679, 0, lcF};
Point(14) = {0.0750,  0.0614, 0, lcF};
Point(15) = {0.0500,  0.0528, 0, lcF};
Point(16) = {0.0250,  0.0407, 0, lcF};
Point(17) = {0.0125,  0.0305, 0, lcF};
Point(18) = {0.0000,  0.0000, 0, lcF};
Point(19) = {0.0125, -0.0205, 0, lcF};
Point(20) = {0.0250, -0.0257, 0, lcF};
Point(21) = {0.0500, -0.0322, 0, lcF};
Point(22) = {0.0750, -0.0360, 0, lcF};
Point(23) = {0.1000, -0.0383, 0, lcF};
Point(24) = {0.1500, -0.0411, 0, lcF};
Point(25) = {0.2000, -0.0422, 0, lcF};
Point(26) = {0.2500, -0.0420, 0, lcF};
Point(27) = {0.3000, -0.0410, 0, lcF};
Point(28) = {0.4000, -0.0376, 0, lcF};
Point(29) = {0.5000, -0.0326, 0, lcF};
Point(30) = {0.6000, -0.0264, 0, lcF};
Point(31) = {0.7000, -0.0196, 0, lcF};
Point(32) = {0.8000, -0.0126, 0, lcF};
Point(33) = {0.9000, -0.0058, 0, lcF};
Point(34) = {0.9500, -0.0025, 0, lcF};

Spline(1) = {1:18};
Spline(2) = {18:34, 1};

Curve Loop(1) = {1, 2};

Point(35) = {-4, -4, 0,      lcC};
Point(36) = {-4,  4, 0,      lcC};
Point(37) = { 5,  4, 0, 10 * lcF};
Point(38) = { 5, -4, 0, 10 * lcF};

Curve(3) = {35, 36};
Curve(4) = {36, 37};
Curve(5) = {37, 38};
Curve(6) = {38, 35};

Curve Loop(2) = {3:6};
Plane Surface(1) = {2, 1};

Field[1] = BoundaryLayer;
Field[1].EdgesList = {1, 2};
Field[1].FanPointsList = {1};
Field[1].hwall_n = lcF / 3;
Field[1].thickness = 0.02;
Field[1].ratio = 1.1;
Field[1].Quads = 1;

BoundaryLayer Field = 1;

Mesh.BoundaryLayerFanElements = 10;

Physical Curve("airfoil") = {1, 2};
Physical Curve("inlet")   = {3};
Physical Curve("top")     = {4};
Physical Curve("outlet")  = {5};
Physical Curve("bottom")  = {6};

Physical Surface("domain") = {1};
