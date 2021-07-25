clear
clc
close all

p = 998; % rho [kgm^-3]
G = 9.81; % acceleration due to gravity [ms^-2]
w = 180.6; % pump rotational speed [rads^-1]
D = 0.060; % pump diameter [m]

V_dot  = [100 200 300 400 500 600 700]';
H = [180 185 175 170 150 95 94]';
eta = [32 54 70 79 81 66 38]';

f = fit(V_dot, H, 'poly2');
g = fit(V_dot, eta, 'poly2');

A = coeffvalues(f);
C = coeffvalues(g);

syms x bhp
g2 = C(1) * x^2 + C(2) * x + C(3); % make function from fitted curve
f2 = A(1) * x^2 + A(2) * x + A(3); % make function from fitted curve

g2dot = diff(g2); % differentiate fitted curve for eta

BEP = solve(g2dot); % solve fitted curve for eta to find max

BEP = double(BEP)

BEP_eta = double(subs(g2, x, BEP))

BEP_head = double(subs(f2, x, BEP))

i = 1;
V_dot_A_array = [100:1:700];
bhp_array = zeros(length(V_dot_A_array), 1);
CP_array = zeros(length(V_dot_A_array), 1);
CQ_array = zeros(length(V_dot_A_array), 1);
CH_array = zeros(length(V_dot_A_array), 1);
eta_array = zeros(length(V_dot_A_array), 1);

g3 = @(x) C(1) * x^2 + C(2) * x + C(3);
f3 = @(x) A(1) * x^2 + A(2) * x + A(3);

eta_A_max = 0;

for V_dot_A = V_dot_A_array
    
eta_A = g3(V_dot_A);
H_A = f3(V_dot_A);

bhp = (p * G * V_dot_A * H_A) / (eta_A/100) * (1/100) ^ 4;

CP = bhp / (p * w^3 * D^5);
CQ = V_dot_A / (w * (D*100)^3);
CH = (G * (H_A/100)) / (w^2 * D^2);

if eta_A> eta_A_max
    eta_A_max = eta_A;
    max_pos = CQ;
end

if V_dot_A == 500
    continue
end

CH_array(i) = CH;
CQ_array(i) = CQ;
CP_array(i) = CP;
bhp_array(i) = bhp;
eta_array(i) = eta_A/100;
i = i + 1;
end

eta_A_max;
max_pos;

figure(1)
yyaxis left
hold on
scatter(V_dot, H, 'xb')
plot(f, '-b')
ylabel('H [cm]')
xlabel('Flow rate [cm^3s^-1]')
yyaxis right
plot(V_dot_A_array, bhp_array, '--r')
plot(g, 'g')
scatter(V_dot, eta, 'g')
xline(BEP, '--k')
ylabel('Pump efficiency [%] or bhp [hp]')
legend('Head (raw data)','Head (fitted curve)', 'bhp',...
    'Efficiency (fitted curve)', 'Efficiency (raw data)',...
    'Best efficiency point')

figure(2)
hold on
xlabel('CQ x 100')
CH_array = CH_array .* 10;
CP_array = CP_array .* 100;
CQ_array = CQ_array .* 100;
max_pos = max_pos * 100;
plot(CQ_array(1:600), CH_array(1:600), 'b')
plot(CQ_array(1:600), eta_array(1:600), '-g')
plot(CQ_array(1:600), CP_array(1:600), '--r')
xline(max_pos, '--k')
legend ('CH x 10',  'Pump efficiency', 'CP x 100', 'Best efficiency point')
