// Gmsh project created on Sun Jun 09 20:11:09 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, -1.0, 0, 1.0};
//+
Point(2) = {5.0, -1.0, 0, 1.0};
//+
Point(3) = {5.0, -1.5, 0, 1.0};
//+
Point(4) = {10.0, -1.5, 0, 1.0};
//+
Point(5) = {10.0, 0, 0, 1.0};
//+
Point(6) = {7.0, 0.0, 0, 1.0};
//+
Point(7) = {5.0, 0.0, 0, 1.0};
//+
Point(8) = {3.0, 0.0, 0, 1.0};
//+
Point(9) = {0.0, 0.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 1};
//+
Curve Loop(1) = {8, 9, 1, 2, 3, 4, 5, 6, 7};
//+
Surface(1) = {1};
//+
Physical Point("bottom") = {1, 2, 3, 4};
//+
Physical Curve("beam") = {7, 6};
//+
Physical Curve("bottom") += {1, 2, 3};
//+
Physical Curve("inlet") = {9};
//+
Physical Curve("damping_in") = {8};
//+
Physical Curve("damping_out") = {5};
//+
Physical Surface("fluid") = {1};