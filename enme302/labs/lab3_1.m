% ENME302 lab 2 question 1
clear
close all
clc

% constants, material properties, values
E = 200 * 10^9;         % elastic modulus, Pa
L1 = 4.5;               % frame 1 length, m
L2 = 4.5;               % frame 2 length, m
L3 = 3;                 % frame 3 length, m
A = 5*10^-4;            % frame area, m^2
I = 1*10^-5;            % frame area moment, m^4
alpha1 = 90;            % angle of frame 1
alpha2 = 0;             % angle of frame 2
alpha3 = -90;           % angle of frame 3

magFactor = 30;
n = 50;


% element 1
K1 = local_frame(E,I,A,L1)
[K1hat, lambda1] = global_frame(K1, alpha1)

% element 2
K2 = local_frame(E,I,A,L2)
[K2hat, lambda2] = global_frame(K2, alpha2)

% element 2
K3 = local_frame(E,I,A,L3)
[K3hat, lambda3] = global_frame(K3, alpha3)

% assembly matrices

A1 = [zeros(3) eye(3); zeros(3, 6)]

A2 = eye(6)

A3 = [zeros(3, 6); eye(3) zeros(3)]
        

% finding K_G
K_G_1 = A1 * K1hat * A1'
K_G_2 = A2 * K2hat * A2'
K_G_3 = A3 * K3hat * A3'

K_G = K_G_1 + K_G_2 + K_G_3
 
 
% equivalent nodal loads
f_eq_1 = addTransverseLVL(313920, L2)

F_eq_1 = lambda2' * f_eq_1

Q_UDL = A2 * F_eq_1


q = K_G\Q_total

D1 = A1' * q
D2 = A2' * q
D3 = A3' * q

d1 = lambda1 * D1
d2 = lambda2 * D2
d3 = lambda3 * D3

f1 = K1 * d1
f2 = K2 * d2
f3 = K3 * d3

F1 = lambda1' * f1
F2 = lambda2' * f2
F3 = lambda3' * f3

plot_deflected_shape(0, 0, 0, L1, d1, magFactor, n, L1, alpha1)
hold on
plot_deflected_shape(0, L2, L1, L1, d2, magFactor, n, L2, alpha2)
plot_deflected_shape(L2, L2, L1, 0, d3, magFactor, n, L3, alpha3)
