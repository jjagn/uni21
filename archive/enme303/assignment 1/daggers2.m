clear
clc
close all

pid1 = import_enme303('jdw144_pid1');

% ------------------------------ CONSTANTS --------------------------------

k_p = 30;
k_d = 30;
k_i = 0.01;

M_c = 1.5;          % mass of the cart (kg)
k_m = 0.0177;       % motor back EMF constant (V/rad/s)
k_g = 3.7;          % gearing ratio
R = 1.5;            % resistance to the motor armature (Ohms)
r = 0.018;          % radius of the pinion (m)
D = 7;              % damping present in the system

B = (k_m * k_g) / (M_c * R * r);
C = D / M_c + (k_m^2 * k_g^2) / (M_c * R * r^2);

num = [B*k_i B*k_d B*k_p];
den = [1+B*k_i C+B*k_d B*k_p];

% ---------------------------- STEP RESPONSE ------------------------------

figure('defaultAxesFontSize', 14)
opts = stepDataOptions('StepAmplitude', 0.1);
[Y] = step(num, den, opts);
t = [0:0.0002611647:1];
plot(t, Y, 'linewidth', 3);

hold on

TF = tf(num, den);
info = stepinfo(TF)
damp(TF)

pid1_x = [0:0.001:5.999];
pid1_y1 = pid1(799:6798, 4);
plot(pid1_x, pid1_y1, '--', 'linewidth', 2);
pid1_y2 = pid1(799:6798, 2);
plot(pid1_x, pid1_y2, 'linewidth', 3);
grid on
legend('Theoretical Position', 'Target Position', 'Actual Position')

ylim([0, 0.18])
xlim([0, 6])

ylabel('Cart Postion (m)')
xlabel('Time (s)')