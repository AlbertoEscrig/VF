// vim: syntax=gmsh

SetFactory("OpenCASCADE");

p1 = newp; Point(p1) = { 0, 0, 0};
p2 = newp; Point(p2) = {20, 0, 0};
p3 = newp; Point(p3) = {20, 2, 0};
p4 = newp; Point(p4) = { 0, 2, 0};

c1 = newc; Curve(c1) = {p1, p2};
c2 = newc; Curve(c2) = {p2, p3};
c3 = newc; Curve(c3) = {p3, p4};
c4 = newc; Curve(c4) = {p4, p1};

c = newcl; Curve Loop(c) = {c1, c2, c3, c4};
s = news; Plane Surface(s) = {c};

Transfinite Curve{c1, c3} = 250; 
Transfinite Curve{c2, c4} = 50 Using Bump 0.3;
Transfinite Surface{:};
Recombine Surface{:};

Physical Curve("walls")  = {c1, c3};
Physical Curve("inlet")  = {c4};
Physical Curve("outlet") = {c2};

Physical Surface("domain") = {s};
