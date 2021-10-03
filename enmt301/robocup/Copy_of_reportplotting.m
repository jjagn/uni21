clear; clc; close all

fsize = 14;
lwidth = 1.5;

steps = [0 2 4 6 8 10 20 30 40 50 60 70 80 90 100]
measured = [0 0 640 610 580 560 550 520 240 200 180 170 165 160 0]

dettime = ones(1,15).* 3

plot(steps, measured, 'LineWidth', lwidth)
hold on; yyaxis right
plot(steps, dettime, 'LineWidth', lwidth)

ylabel("Detection time [ms]", "FontSize", fsize)

yyaxis left
ylabel("Sensor value", "FontSize", fsize)
xlabel("Distance [cm]", "FontSize", fsize)

legend("Sensor value", "Detection time", "FontSize", fsize)