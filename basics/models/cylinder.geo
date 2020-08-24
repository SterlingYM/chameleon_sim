Merge "cylinder.step";
Physical Surface("inner", 1) = {4};
//+
Physical Surface("outer", 2) = {1, 2, 3};
//+
Physical Volume("vacuum", 3) = {1};
//+
Transfinite Line {4} = 30 Using Progression 1;
//+
Transfinite Line {2, 1, 3} = 50 Using Progression 1;
