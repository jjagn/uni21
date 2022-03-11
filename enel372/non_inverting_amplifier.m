%%% NON-INVERTING OPAMP NOISE FACTOR
clear; close all; clc

%% DECLARING VARIABLES
R1 = 1E3;
R2 = 4E3;
RS = 10E3;
in = 0.01E-12;
ia = in;
en = 7E-9;
fkt = 1.62E-20;

%% CALCULATIONS
RP = R1*R2/(R1+R2)

ea2 = en^2+fkt*RP+(in*RP)^2

NF = 10*log10(1+(ea2+(ia*RS)^2)/fkt/RS)

