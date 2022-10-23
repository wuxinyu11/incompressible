
L = 1.0;
n = 96;

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};

Line(1) = {1,2};

Transfinite Curve{1} = n+1;
Physical Curve("Ω") = {1};
Physical Point("Γᵍ") = {1,2};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 1;