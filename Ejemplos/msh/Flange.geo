// vim: syntax=gmsh

SetFactory("Built-in");

Merge "flange.stl";

Surface Loop(1) = Surface{:};
Volume(1) = {1};

Physical Volume("domain") = {1};

Physical Surface("hot")  = {4};
Physical Surface("cold") = {2};
Physical Surface("rest") = {1, 3};
