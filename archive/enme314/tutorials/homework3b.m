% solving system of equations for enme314 homework 3 (DAN ZHAO HW 1)
clear
close all
clc

syms hl v1 v2 v3 f1 f2 f3 vdot vdot1 vdot2 vdot3 hl1 hl2 hl3 hpump re1 re2 re3

p = 999.1;
mu = 1.138 * 10 ^ -3;
g = 9.81;
Welec = 8000;
D1 = 0.03;
D2 = 0.05;
D3 = 0.04;
L = 25.1; % 25.1
eta = 0.68;
eps = 0.3/1000; % 0
zA = 2.1;
zB = 9.1; % 9.1

eqn1 = Welec == p*vdot*g*hpump/eta;
eqn2 = hl1 == hl2;
eqn3 = hl2 == hl3;
eqn4 = hl == hl1;
eqn5 = hpump == (zB-zA)+hl;
eqn6 = v1 == vdot1/(pi*D1^2/4);
eqn7 = v2 == vdot2/(pi*D2^2/4);
eqn8 = v3 == vdot3/(pi*D3^2/4);
eqn9 = vdot == vdot1 + vdot2 + vdot3;
eqn10 = re1 == p * v1 * D1 / mu;
eqn11 = re2 == p * v2 * D2 / mu;
eqn12 = re3 == p * v3 * D3 / mu;
eqn13 = 1/sqrt(f1) == -2.0 * log10((eps/D1)/3.7 + 2.51/(re1 * sqrt(f1)));
eqn14 = 1/sqrt(f2) == -2.0 * log10((eps/D2)/3.7 + 2.51/((re2 * sqrt(f2))));
eqn15 = 1/sqrt(f3) == -2.0 * log10((eps/D3)/3.7 + 2.51/((re3 * sqrt(f3))));
eqn16 = hl1 == f1*L/D1*(v1^2)/(2*g);
eqn17 = hl2 == f2*L/D2*(v2^2)/(2*g);
eqn18 = hl3 == f3*L/D3*(v3^2)/(2*g);

eqns = [eqn1; eqn2; eqn3; eqn4; eqn5; eqn6; eqn7; eqn8; eqn9; eqn10; ...
    eqn11; eqn12; eqn13; eqn14; eqn15; eqn16; eqn17; eqn18];

vars = [hl v1 v2 v3 f1 f2 f3 vdot vdot1 vdot2 vdot3 hl1 hl2 hl3 hpump...
    re1 re2 re3];

S = vpasolve(eqns, vars);

hl = double(S.hl)
v1 = double(S.v1)
f1 = double(S.f1)
f2 = double(S.f2)
f3 = double(S.f3)
vdot = double(S.vdot)
vdot1 = double(S.vdot1)
vdot2 = double(S.vdot2)
vdot3 = double(S.vdot3)
hl1 = double(S.hl1)
hl2 = double(S.hl2)
hl3 = double(S.hl3)
hpump = double(S.hpump)
re1 = double(S.re1)
re2 = double(S.re2)
re3 = double(S.re3)