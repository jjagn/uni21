clear

close all

clc

Cd = 0.98;
A1 = 338.6*10^-6;
A2 = 84.6*10^-6;

rho = 998.2;

a = mean([21.9 21.4 20.7]);
b = mean([15.3 14.8 14.5]);
c = mean([9.7 9.6 9.3]);
d = mean([6.0 5.6 5.4]);
e = mean([2.5 2.2 2.3]);
f = mean([0.53 0.45 0.44]);

pressureDrop = [a b c d e f];

pressureDrop = pressureDrop*100

pressureDrop = pressureDrop'

Q = Cd * A2 * sqrt((2.*pressureDrop) ./ (rho*(1-A2^2/A1^2)))