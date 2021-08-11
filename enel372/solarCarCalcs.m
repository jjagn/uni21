clear
close all
clc

Vs = 17.5; % may be 14.5
Vo = 2;
Is = 0.26;
Io = 2;
% fs = 10 * 10^3;
fs = linspace(10*10^3, 100*10^3, 500);
% Ipp = ;
% Vpp = ;
% D = linspace(0.36, 0.83, 100);
D = .5;

L = (Vs .* D .* (1 - D)) ./ (2 .* fs .* Io);
    
plot(fs, L)