clear
close all
clc

Vs = 15; % may be 14.5
Vo = 7;
Is = 0.6;
Io = 2.5;
fs = 50 * 10^3;
% fs = linspace(10*10^3, 100*10^3, 500);
% Ipp = ;
% Vpp = ;
% D = linspace(0.36, 0.83, 100);
D = .5;

Lmin = Vs / (8 * fs * Io)
dV = Vo/10;
C = Io/(4*fs*dV) * 1e03

Lmin * 1e03


Ipp = Vs * (1-D) / (fs - L)