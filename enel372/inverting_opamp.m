%% INVERTING OPAMP NOISE FACTOR
clear; close all; clc

%% DECLARING VARIABLES
R1 = 20E3;
R2 = 100E3;
RS = 20E3;
in = 0.06E-12;
en = 7E-9;
fkt = 1.62E-20;

%% CALCULATIONS
ia2 = in^2 + fkt * 1/R2 + fkt * 1/R1

ea2 = en^2+R1^2*ia2

NF = 10*log10(1+(ea2+(ia2*RS^2))/fkt/RS)

