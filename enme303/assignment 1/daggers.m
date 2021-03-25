clear
clc
close all

% ------------------------------ CONSTANTS --------------------------------

k_d = 30;
k_p = 30;
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
den = [1+B*k_i C+B*k_d B*k_d];

% ---------------------------- STEP RESPONSE ------------------------------

opts = stepDataOptions('StepAmplitude', 0.1);
[Y] = step(num, den, opts);
step(num, den, opts);
sys = tf(num, den);

info = stepinfo(sys)

% ------------------------------- VOLTAGE ---------------------------------

E = 0.1 - Y;
E_dot = [diff(E); 0];
E_d_dot = trapz(E);

V = E*k_p + E_dot*k_d + E_d_dot*k_i;

figure(2)
plot(V)