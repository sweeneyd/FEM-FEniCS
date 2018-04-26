// Gmsh project created on Mon Oct 23 12:29:11 2017
SetFactory("OpenCASCADE");
edgeMeshSize = 0.1;
//+
Point(1) = {0, 0, 0, edgeMeshSize};
//+
Point(2) = {1, 0, 0, edgeMeshSize};
//+
Point(3) = {1, 1, 0, edgeMeshSize};
//+
Point(4) = {0, 1, 0, edgeMeshSize};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {0.5, 0.5, 0, 0.1, 0, 2*Pi};
//+ Add this to change mesh size
Characteristic Length {5} = edgeMeshSize/10;
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Line Loop(2) = {5};
//+
Plane Surface(1) = {2, 1};
//+
Physical Line(1) = {1};
//+
Physical Line(2) = {2};
//+
Physical Line(3) = {3};
//+
Physical Line(4) = {4};
//+
Physical Line(5) = {5};
//+
Physical Surface(0) = {1};
