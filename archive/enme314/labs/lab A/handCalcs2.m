clear
close all
clc

p_ambient = 102138;
t_ambient_C = 22;
t_ambient_K = t_ambient_C + 273.15;

Ryd = 287; % Rydberg constant for ideal gas equation

g = 9.81;
x3 = 0.35; % metres

F = [0 0 0 0 0 0 0 0 0 -1 -1 -1]; % 0 = flat, 1 = hemisphere

m = [151.9 151.9 151.9 151.9 151.9 151.9 152 152 152 304 304 304];
m = m/1000;

x2 = [259 223 189 286 291 316 194 218 216 174 161 151];
x2 = x2/1000;

h = [100 180 260 100 180 260 270 180 100 100 180 270];

d = [70 70 70 150 150 150 110 110 110 110 110 110];

jet_temp_C = [34 34 34 34.3 34.4 34.5 36.7 36.7 37 36.3 36.6 37.1];

pressure_bell_mmH20 = [58.3 58.3 58.2 58.2 58.2 58.2 65 63.9 63.9 65.3 65.1 65.1];

jet_temp_K = jet_temp_C + 273.15;
pressure_bell_Pa = pressure_bell_mmH20 * 9.80665;


data = [pressure_bell_Pa', jet_temp_K', x2', m', F'];

A2 = (34/1000)^2*pi/4;

A3 = (39/1000)^2*pi/4;

Cd = 0.98;

rho1 = p_ambient/(Ryd*t_ambient_K);

rho2 = (p_ambient - pressure_bell_Pa) ./ (Ryd .* jet_temp_K);

dp = p_ambient' - pressure_bell_Pa'

rho2 = rho2';

m_dot = Cd .* rho2 .* A2 .* sqrt((2.*pressure_bell_Pa') ./rho1);

v2 = sqrt(2*pressure_bell_Pa'./rho2)
% v2b = m_dot ./ (rho2 .* A2)

v3 = A2/A3 .* v2;
% v3b = A2/A3 .* v2b;

v4 = v3 .* F';
% v4b = v3b .* F';

F_D_E = m' .* g .* (x2'./x3');

F_D_T = m_dot .* (v3 - v4);
% F_D_T_b = m_dot .* (v3b - v4b)

data2 = [data(:, 1:4) F_D_E F_D_T];

p_unc_abs = 0.1 * 9.80665;
p1_unc_abs = 1;
t_unc_abs = 0.2;
t1_unc_abs = 0.2;
x2_unc_abs = 3*10^-3;
x3_unc_abs = 1*10^-3;
m_unc_abs = 5*10^-5;
d2_unc_abs = 0.5;
d3_unc_abs = 0.5;

p_unc_rel = p_unc_abs ./ pressure_bell_Pa';
p1_unc_rel = p1_unc_abs / p_ambient;
t_unc_rel = t_unc_abs ./ jet_temp_C';
t1_unc_rel = t_unc_abs / 22;
m_unc_rel = m_unc_abs ./ m';
x2_unc_rel = x2_unc_abs ./ x2';
x3_unc_rel = x3_unc_abs / x3;
d2_unc_rel = d2_unc_abs / 34;
d3_unc_rel = d3_unc_abs / 39;

relative_uncertainty_experimental = x2_unc_rel + x3_unc_rel + m_unc_rel;
absolute_uncertainty_experimental = F_D_E .* relative_uncertainty_experimental;

relative_uncertainty_theoretical = 2.5 * p_unc_rel + 1.5 * t_unc_rel + 1.5 * p1_unc_rel + 4 * d2_unc_rel + 4 * d3_unc_rel + 0.5 * t_unc_rel;
absolute_uncertainty_theoretical = relative_uncertainty_theoretical .* F_D_T;

Experimental_Drag_data = [F_D_E, absolute_uncertainty_experimental]

Theoretical_Drag_data = [F_D_T, absolute_uncertainty_theoretical]


data3 = [h' d' Experimental_Drag_data Theoretical_Drag_data]

