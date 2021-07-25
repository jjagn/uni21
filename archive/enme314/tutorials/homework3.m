% solving system of equations for enme314 homework 3 (DAN ZHAO HW 1)
clear
close all
clc

syms hl v1 f1 f2 f3

p = 999.1;
mu = 1.138 * 10 ^ -3;
g = 9.81;
Welec = 8000;
D1 = 0.03;
D2 = 0.05;
D3 = 0.04;
L = 25.1; % 25.1
eta = 0.68;
eps = 0.3 / 1000; % 0
zA = 2.1;
zB = 9.1; % 9.1

v2 = sqrt(f1/f2 * D2/D1) * v1; 
v3 = sqrt(f1/f3 * D3/D1) * v1; 

eqn1 = Welec  == (p * g * pi/4 * (D1 ^ 2 * v1 + D2 ^ 2 * v2 + D3 ^ 2 * v3) * ((zB - zA) + hl))/eta;
eqn2 = hl == f1 * L/D1 * v1^2/(2*g);
eqn3 = 1/sqrt(f1) == -2.0 * log10((eps/D1)/3.7 + 2.51/((p * v1 * D1 * sqrt(f1))/mu));
eqn4 = 1/sqrt(f2) == -2.0 * log10((eps/D2)/3.7 + 2.51/((p * v2 * D2 * sqrt(f2))/mu));
eqn5 = 1/sqrt(f3) == -2.0 * log10((eps/D3)/3.7 + 2.51/((p * v3 * D3 * sqrt(f3))/mu));

eqns = [eqn1; eqn2; eqn3; eqn4; eqn5];

vars = [hl, v1, f1, f2, f3];

initial = [18 5 0.015 0.015 0.015];

S = vpasolve(eqns, vars, initial);

hl = double(S.hl);
v1 = double(S.v1);
f1 = double(S.f1);
f2 = double(S.f2);
f3 = double(S.f3);

v2 = sqrt(f1/f2 * D2/D1) * double(S.v1); 
v3 = sqrt(f1/f3 * D3/D1) * double(S.v1);

v_dot1 = pi*D1^2/4*double(S.v1);
v_dot2 = pi*D2^2/4*double(v2);
v_dot3 = pi*D3^2/4*double(v3);

v_dot = v_dot1 + v_dot2 + v_dot3;

ans = [v_dot v_dot1 v_dot2 v_dot3 hl]