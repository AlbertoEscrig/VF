// vim: syntax=gmsh

SetFactory("OpenCASCADE");

mm = 1e-3;

P1 = newp; Point(P1) = {206 * mm, 0, 0};
P2 = newp; Point(P2) = {290 * mm, 0, 0};
P3 = newp; Point(P3) = {290 * mm, -33.25 * mm / 2, 0};
P4 = newp; Point(P4) = {206 * mm, -50.8 * mm / 2, 0};

L1 = newc; Curve(L1) = {P1, P2};
L2 = newc; Curve(L2) = {P2, P3};
L3 = newc; Curve(L3) = {P3, P4};
L4 = newc; Curve(L4) = {P4, P1};

LL1 = newcl; Curve Loop(LL1) = {L1, L2, L3, L4};

S1 = news; Plane Surface(S1) = {LL1};

P5 = newp; Point(P5) = {290 * mm, 33.25 * mm / 2, 0};
P6 = newp; Point(P6) = {206 * mm, 50.8 * mm / 2, 0};

L5 = newc; Curve(L5) = {P2, P5};
L6 = newc; Curve(L6) = {P5, P6};
L7 = newc; Curve(L7) = {P6, P1};

LL2 = newcl; Curve Loop(LL2) = {L1, L5, L6, L7};

S2 = news; Plane Surface(S2) = {LL2};

Transfinite Curve{L1, L3, L6} = 75;
Transfinite Curve{L2, L4, L5, L7} = 15;

Transfinite Surface{S1, S2};
Recombine Surface{S1, S2};

out[] = Extrude {-206 * mm, 0, 0} { Curve{L4, L7}; Layers{150}; Recombine; };

L8  = out[0];
L9  = out[2];
L10 = out[4];
L11 = out[5];

out[] = Extrude {-20.6 * mm, 0, 0} { Curve{L10}; Layers{15}; Recombine; };

L12 = out[0];
L13 = out[2];
L14 = out[3];

Physical Curve("inlet")    = {L12};
Physical Curve("outlet")   = {L2, L5};
Physical Curve("wall")     = {L3, L6, L8, L9, L11, L13, L14};
Physical Surface("domain") = Surface{:};
