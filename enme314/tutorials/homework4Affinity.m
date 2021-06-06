clear
clc
close all

pA = 998; % rho for water [kgm^-3]
pB = 1226; % rho for R-134a [kgm^-3]
G = 9.81; % acceleration due to gravity [ms^-2]
wA = 180.6; % pump rotational speed [rads^-1]
dA = 0.060; % pump diameter [m]

% required values for pump B
BEP_B = 2400; % best efficiency point for pump B, [cm^3s^-1]
BEP_head_B = 450; % available head from pump B at BEP [cm]


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

BEP_A = solve(g2dot); % solve fitted curve for eta to find max

BEP_A = double(BEP_A)

BEP_eta = double(subs(g2, x, BEP_A))

BEP_head_A = double(subs(f2, x, BEP_A))

V_dot_A_array = [100:1:700];

g3 = @(x) C(1) * x^2 + C(2) * x + C(3);
f3 = @(x) A(1) * x^2 + A(2) * x + A(3);

eta_A_max = 0;

for V_dot_A = V_dot_A_array
    
eta_A = g3(V_dot_A);
H_A = f3(V_dot_A);

if eta_A > eta_A_max
    eta_A_max = eta_A;
    bhp_A = (pA * G * V_dot_A * H_A) / eta_A * (1/100) ^ 3;
end


end

syms bhp_B wB dB

eqn1 = BEP_B == BEP_A * (wB/wA) * (dA/dB)^3;
eqn2 = BEP_head_B == BEP_head_A * (wB/wA)^2 * (dB/dA)^2;
eqn3 = bhp_B == bhp_A * pB/pA * (wB/wA)^3 * (dB/dA)^5;

eqn4 = BEP_B/BEP_A == (wB/wA) * (dB/dA)^3;
eqn5 = BEP_head_B/BEP_head_A == (wB/wA)^2 * (dB/dA)^2;
eqn6 = bhp_B/bhp_A == pB/pA * (wB/wA)^3 * (dB/dA)^5;


eqns = [eqn1, eqn2, eqn3];
vars = [bhp_B, dB, wB];

eqns2 = [eqn4 eqn5 eqn6];

I = [0 0 0];

S = vpasolve(eqns, vars, I);
S2 = vpasolve(eqns2, vars, I);

bhp_B_soln = (double(S.bhp_B));
wB_soln = (double(S.wB));
dB_soln = (double(S.dB));

%solns = [bhp_B_soln wB_soln dB_soln];
solns = real([bhp_B_soln wB_soln dB_soln])

bhp_B_soln2 = (double(S2.bhp_B));
wB_soln2 = (double(S2.wB));
dB_soln2 = (double(S2.dB));

%solns3 = [bhp_B_soln2 wB_soln2 dB_soln2];
solns2 = real([bhp_B_soln2 wB_soln2 dB_soln2])

