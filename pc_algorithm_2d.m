function [sign,Z_positive,Z_negative] = pc_algorithm_2d(NUM_CNST, NUM_DV, loop, xval, xmin, xmax, xold1, xold2, low, upp, ak, a, c, d)

m0 = -1;
n0 = -1;
a0 = 0;
b0 = 0;

height = 0.2;
steep = -100;

m1 = height;
n1 = steep;
a1 = -0.7;
b1 = -0.1;

m2 = height;
n2 = steep;
a2 = -0.7;
b2 = 0.1;

m3 = height;
n3 = steep;
a3 = -0.5;
b3 = -0.1;

m4 = height;
n4 = steep;
a4 = -0.5;
b4 = 0.1;

objF = @(x) m0*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
    + m1*exp(n1*((x(1)-a1)^2+(x(2)-b1)^2))...
    + m2*exp(n2*((x(1)-a2)^2+(x(2)-b2)^2))...
    + m3*exp(n3*((x(1)-a3)^2+(x(2)-b3)^2))...
    + m4*exp(n4*((x(1)-a4)^2+(x(2)-b4)^2));

gradF = @(x) [2*m0*n0*(x(1)-a0)*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
    + 2*m1*n1*(x(1)-a1)*exp(n1*((x(1)-a1)^2+(x(2)-b1)^2))...
    + 2*m2*n2*(x(1)-a2)*exp(n2*((x(1)-a2)^2+(x(2)-b2)^2))...
    + 2*m3*n3*(x(1)-a3)*exp(n3*((x(1)-a3)^2+(x(2)-b3)^2))...
    + 2*m4*n4*(x(1)-a4)*exp(n4*((x(1)-a4)^2+(x(2)-b4)^2));
    2*m0*n0*(x(2)-b0)*exp(n0*((x(1)-a0)^2+(x(2)-b0)^2))...
    + 2*m1*n1*(x(2)-b1)*exp(n1*((x(1)-a1)^2+(x(2)-b1)^2))...
    + 2*m2*n2*(x(2)-b2)*exp(n2*((x(1)-a2)^2+(x(2)-b2)^2))...
    + 2*m3*n3*(x(2)-b3)*exp(n3*((x(1)-a3)^2+(x(2)-b3)^2))...
    + 2*m4*n4*(x(2)-b4)*exp(n4*((x(1)-a4)^2+(x(2)-b4)^2))];

x_cur = xval;

g_cur = objF(x_cur);

h_positive = gradF(x_cur);
h_negative = -gradF(x_cur);


% Obtainting Z+
obj = g_cur;
Dobj = h_positive;
f0val = obj; 
df0dx = Dobj; 
df0dx2 = 0;
cnst01 = objF(x_cur) -500;
Dcnst01 = gradF(x_cur)';
fval = cnst01';
dfdx = Dcnst01;
dfdx2 = 0;

[xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp] = ...
        mmasub(NUM_CNST, NUM_DV, loop, x_cur, xmin, xmax, xold1, xold2, f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, ak, a, c, d);

x_new_positive = xmma;
g_new_positive = objF(x_new_positive);
Z_positive = g_new_positive - g_cur;


% Obtainting Z-
obj = g_cur;
Dobj = h_negative;
f0val = obj; 
df0dx = Dobj; 
df0dx2 = 0;
cnst01 = objF(x_cur) -500;
Dcnst01 = gradF(x_cur)';
fval = cnst01';
dfdx = Dcnst01;
dfdx2 = 0;

[xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp] = ...
        mmasub(NUM_CNST, NUM_DV, loop, x_cur, xmin, xmax, xold1, xold2, f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, ak, a, c, d);

x_new_negative = xmma;
g_new_negative = objF(x_new_negative);
Z_negative = g_new_negative - g_cur;


% Discriminants
if (Z_positive - Z_negative < 0)
    sign = 1;
else
    sign = -1;
end





