a = 1.0;
b = 1.0;
n = 74;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  a,   b, 0.0};
Point(4) = {0.0,   b, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(5) = {1,2,3,4};

Plane Surface(1) = {5};
Transfinite Curve{1,2,3,4} = n+1;
//Transfinite Curve{1,3} = n;
//Transfinite Curve{2,4} = n+1;
Transfinite Surface{1};

Physical Curve("Γᵍ") = {1,2,3,4};
Physical Curve("Γᵗ") = {};
Physical Surface("Ω") = {1};

//Mesh.ElementOrder = 2;
Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;