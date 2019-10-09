/* Objective function */
min: x1c + x4c;

x1 + x3 = 1;
x2 + x4 = 1;
-0.9 x3 + 0.1 x2 = 0;
-0.1 x2 + 0.9 x3 = 0;
0.9 x1 + 0.1 x2 = 0.9;
0.9 x3 + 0.1 x4 = 0.1;
x1 - x1c <= 0.5;
x1 + x1c >= 0.5;
x4 - x4c <= 0.5;
x4 + x4c >= 0.5;


/* Variable bounds */