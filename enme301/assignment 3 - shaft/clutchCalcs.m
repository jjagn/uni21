clear
close all
clc

T0 = 17362; % torque transmitted through clutch
Tw = T0; Tp = T0;

D = 0.5;
D_d_array = (0.01:0.01:1); % starting at 50mm because anything below that makes actuation force crazy
F_array_Tw = zeros(5, length(D_d_array));
F_array_Tp = zeros(5, length(D_d_array));
i = 1;

while i <= length(D_d_array)
 % let D = maximum allowable 500mm
d = D_d_array(i) * D; % let d = 0.5D
f = 0.55; % let friction coefficient be that for standard clean and dry steel on steel
for n = (1:5)
F_Tw = ((Tw * 4) / ((D + d) * f)) / n;
F_Tp = ((Tp * 3) / f * (D^2 - d^2)/(D^3 - d^3)) / n;

F_array_Tw(n, i) = F_Tw;
F_array_Tp(n, i) = F_Tp;
end
i = i + 1;
end

figure(1)
plot(D_d_array, F_array_Tw, 'LineWidth', 2)
title('Constant wear')
xlabel('d/D')
ylabel('Actuation force [N]')
legend('1 pair of plates', '2 pairs', '3 pairs', '4 pairs', '5 pairs')

figure(2)
plot(D_d_array, F_array_Tp, 'LineWidth', 2)
title('Constant pressure')
xlabel('d/D')
ylabel('Actuation force [N]')
legend('1 pair of plates', '2 pairs', '3 pairs', '4 pairs', '5 pairs')

figure(3)
plot(D_d_array, F_array_Tw(1, :), 'LineWidth', 1)
hold on
plot(D_d_array, F_array_Tp(1, :), 'LineWidth', 1)
title('Constant wear vs constant pressure, 1 plate')
xlabel('d/D')
ylabel('Actuation force [N]')
legend('Constant wear', 'Constant pressure')

