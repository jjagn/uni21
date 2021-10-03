clear; clc; close all

fsize = 14;
lwidth = 1.5;

steps = linspace(100,1000,10)
calc = linspace(2,20,10)
measured = [1.9 3.9 5.9 7.8 9.8 11.7 13.6 15.6 17.5 19.5]


plot(steps, calc, 'LineWidth', lwidth)
hold on
plot(steps, measured, 'LineWidth', lwidth)

xlabel("Number of steps", "FontSize", fsize)
ylabel("Distance traveled [mm]", "FontSize", fsize)

legend("Calculated distance", "Measured distance", "FontSize", fsize)