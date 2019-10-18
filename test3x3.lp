/* Objective function */
min: x1c + x2c + x3c + x4c + x5c + x6c + x7c + x8c + x9c;

x1 - x1c <= 0.5;
x1 + x1c >= 0.5;

x2 - x2c <= 1;
x2 + x2c >= 1;

x3 - x3c <= 1/3;
x3 + x3c >= 1/3;

x4c = 0;
x5c = 0;

x6 - x6c <= 1/3.;
x6 + x6c >= 1/3.;

x7 - x7c <= 0.5;
x7 + x7c >= 0.5;

x8c = 0;

x7 - x7c <= 1/3.;
x7 + x7c >= 1/3.;

x1 + x4 + x7 = 1;
x2 + x5 + x8 = 1;
x3 + x6 + x9 = 1;

-0.3 x7 - 0.3 x4 + 0.3 x2 + 0.4 x3 = 0;

0.3 x1 + 0.3 x2 + 0.4 * x3 = 0.3;
0.3 x4 + 0.3 x5 + 0.4 * x6 = 0.3;
0.3 x7 + 0.3 x8 + 0.4 * x9 = 0.4;


/* Variable bounds */