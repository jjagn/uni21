% ENME302 lab 2 question 1
clear
close all
clc

% format long

% constants, material properties, values
E = 200 * 10^9;         % elastic modulus, Pa
L1 = 5;                  % frame length, m
L2 = 4;

A = 6*10^-4;            % frame area, m^2
I = 8*10^-6;            % frame area moment, m^4
alpha1 = 53.1;             % angle of frame 1
alpha2 = -90;             % angle of frame 2

magFactor = 30;
n = 50;


% element 1
K1 = local_frame(E,I,A,L1)
[K1hat, lambda1] = global_frame(K1, alpha1)

% element 2
K2 = local_frame(E,I,A,L2)
[K2hat, lambda2] = global_frame(K2, alpha2)

% assembly matrices
A1 = [0 0 0 1 0 0;...
      0 0 0 0 1 0;...
      0 0 0 0 0 1]
  
A2 = [1 0 0 0 0 0;...
      0 1 0 0 0 0;...
      0 0 1 0 0 0]
  

% finding K_G
K_G_1 = A1 * K1hat * A1'
K_G_2 = A2 * K2hat * A2'

K_G = K_G_1 + K_G_2

% external reactions
f_eq_1 = addTransverseLVL(-78480, L2)

F_eq_1 = lambda2' * f_eq_1

Q_LVL = A2 * F_eq_1

format long

q = K_G\Q_LVL

format short

D1 = A1' * q
D2 = A2' * q

d1 = lambda1 * D1
d2 = lambda2 * D2

f1 = K1 * d1
f2 = K2 * d2

F1 = lambda1' * f1
F2 = lambda2' * f2;
F2 - F_eq_1

eps1 = (d1(4) - d1(1))/L1
eps2 = (d2(4) - d2(1))/(L2)

plot_deflected_shape(0, 0, 3, 4, d1, magFactor, n, L1, 53.1)
hold on
plot_deflected_shape(3, 4, 3, 0, d2, magFactor, n, L2, -90)
axis equal
