
SetFactory("OpenCASCADE");

mm = 1e-3;

P1 = newp; Point(P1) = {206 * mm, 0,0};
P2 = newp; Point(P2) = {290 * mm, 0,0};
P3 = newp; Point(P3) = {290 * mm, -33.25 * mm / 2, 0};
P4 = newp; Point(P4) = {206 * mm, -50.8 * mm / 2, 0};

L1 = newl; Line(L1) = {P1, P2};
L2 = newl; Line(L2) = {P2, P3};
L3 = newl; Line(L3) = {P3, P4};
L4 = newl; Line(L4) = {P4, P1};

LL1 = newll; Line Loop(LL1) = {L1, L2, L3, L4};

S1 = news; Plane Surface(S1) = {LL1};

P5 = newp; Point(P5) = {290 * mm, 33.25 * mm / 2, 0};
P6 = newp; Point(P6) = {206 * mm, 50.8 * mm / 2, 0};

L5 = newl; Line(L5) = {P2, P5};
L6 = newl; Line(L6) = {P5, P6};
L7 = newl; Line(L7) = {P6, P1};

LL2 = newll; Line Loop(LL2) = {L1, L5, L6, L7};

S2 = news; Plane Surface(S2) = {LL2};

Transfinite Line{L1, L3, L6} = 50;
Transfinite Line{L2, L4, L5, L7} = 15;

Transfinite Surface{S1, S2};
Recombine Surface{S1, S2};

Extrude {-206 * mm, 0, 0} { Line{L4, L7}; Layers{100}; Recombine; }

Extrude {-20.6 * mm, 0, 0} { Line{14}; Layers{15}; Recombine; }

Physical Line("inlet") = {17};
Physical Line("outlet") = {7, 2};
Physical Line("wall") = {15, 13, 8, 16, 12, 10, 3};
Physical Surface("dominio") = {14, 13, 12, 11, 6};
