% jackson crawford 
% enme303 assignment 1

clear
clc
close all

pd1 = import_enme303('jdw144_pd1');

M = 1.5;                % mass of cart, kg
km = 0.017;             % motor back emf constant, V/rad/s
kg = 3.7;               % gearing ratio
R = 1.5;                % resistance of motor armature, ohms
r = 0.018;              % radius of pinion, m
D = 7;                  % damping present in physical system, e.g friction

% constants so I don't have to type all that junk out
B = (km * kg) / (M * R * r);
C = D / M + (km^2 * kg^2) / (M * R * r^2);

% controller gains
Kd = 15;             % 15.36
Kp = 120;

% fraction declaration
num = [B*Kd B*Kp];
den = [1 C+B*Kd B*Kp];

system = tf(num, den);

[Wn, Z, P] = damp(system);

opts = stepDataOptions('StepAmplitude', 0.1);

[Y, T] = step(system, opts);

% step(num, den, opts)
% figure()
% step(num, den, opts)
% grid ON

goalLine = ones(length(T)) * 0.1;


figure('DefaultAxesFontSize', 14) 
plot(T, Y, 'LineWidth', 1);
hold on
plot(T, goalLine, ':', 'LineWidth', 1)
xlabel('Time (s)')
ylabel('Amplitude (m)')
xlim([0, 1])
grid ON

% figure()
% rlocus(system)

% % error terms
% E = 0.1-Y;
% E_dot = diff(E);
% E_dot = [0; E_dot];
% 
% % calculating requested voltage using error
% V = Kp*E + Kd*E_dot;
% 
% % scale_factor = max(V)/max(Y);
% 
% % V_scaled = V/scale_factor + 1;
% V_scaled = V + 1;
% 
% figure();
% hold on
% title('System response and voltage requested');
% yyaxis right
% plot(Y)
% ylabel('system response');
% yyaxis left
% plot(V)
% hold off
% 
% % titlestr1 = sprintf('Requested voltage (scaled) by factor of %.3g (3 s.f)', scale_factor);
% titlestr1 = 'Requested voltage';
% ylabel(titlestr1);

info = stepinfo(system);

% 784

figure('DefaultAxesFontSize', 14)
% step(num, den, opts);
hold on
plot(pd1(:, 1), pd1(:, [2,4]), 'LineWidth', 1);
grid ON
title('Cart position and commanded position, P = 120, D = 15.36')
xlim([0, 10]);
legend('cart position', 'position command')
