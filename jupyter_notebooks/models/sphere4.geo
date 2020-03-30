Merge "sphere2.step";
Physical Surface("outer", 1) = {1};
//+
Physical Volume("vacuum", 4) = {1};
//+
Point(3) = {0, 0, 0, 0.01};
//+
Transfinite Line {1} = 10 Using Progression 1;
//+
Transfinite Line {1} = 30 Using Progression 1;
//+
SetFactory("OpenCASCADE");
Sphere(2) = {0, 0, 0, 10, -Pi/2, Pi/2, 2*Pi};
//+
Transfinite Surface {2};
