clear

close all

clc

Cd = 0.7;

D1 = 28.5/1000;
D2 = 18.5/1000;

A1 = (D1^2*pi)/4;
A2 = (D2^2*pi)/4;

rho = 998.2;

a = mean([3.5 3.42 3.39]);
b = mean([2.32 2.25 2.24]);
c = mean([1.34 1.32 1.27]);
d = mean([0.57 0.56 0.55]);
e = mean([0.25 0.25 0.24]);
f = mean([0.1 0.1 0.13]);

pressureDrop = [a b c d e f];

F = [600 500 400 300 200 100]' * (1000*3600)^-1

pressureDrop = pressureDrop*100

pressureDrop = pressureDrop'

Q = Cd * A2 * sqrt((2.*pressureDrop) ./ (rho*(1-A2^2/A1^2)))