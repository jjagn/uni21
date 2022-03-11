%% SALLEN_KEY HPF
clear; close all; clc

syms w0 R3 R4 C1 C2 Q

eqn1 = w0/Q == (C1+C2)/(R4*C1*C2);
eqn2 = w0^2 == 1/(R3*R4*C1*C2);

C1_val = 0.6E-6;
C2_val = C1_val;
cutoff_frequency = 2E3;
w0_val = 2*pi*cutoff_frequency;
Q_val = 1/sqrt(2);
%% SOLVING
eqn3 = solve(eqn1, R4)

R4_val = double(subs(eqn3, {C1,C2,w0,Q}, [C1_val C2_val w0_val Q_val]))


eqn4 = solve(eqn2,R3)
R3 = subs(eqn4, {C1,C2,w0,R4}, [C1_val C2_val w0_val,R4_val]);
R3_val = double(R3)