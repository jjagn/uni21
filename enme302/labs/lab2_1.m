% ENME302 lab 2 question 1
clear
close all
clc

% constants, material properties, values
E = 200 * 10^9;         % elastic modulus, Pa
L = 10;                 % frame length, m
Q1 = 0;                 % reactions, N
Q2 = 140 * 10^3;        % "" ""
Q3 = 0;                 % "" ""
A = 1*10^-5;            % frame area, m^2
I = 5*10^-4;            % frame area moment, m^4
alpha1 = -90;             % angle of frame 1
alpha2 = 0;             % angle of frame 2

magFactor = 100;
n = 50;


% element 1
K1 = local_frame(E,I,A,L)
[K1hat, lambda1] = global_frame(K1, alpha1)

% element 2
K2 = local_frame(E,I,A,L)
[K2hat, lambda2] = global_frame(K2, alpha2)

% assembly matrices
A1 = [0 0 0 1 0 0;...
      0 0 0 0 0 1;...
      0 0 0 0 0 0];
  
A2 = [1 0 0 0 0 0;...
      0 0 1 0 0 0;...
      0 0 0 0 0 1];

% finding K_G
K_G_1 = A1 * K1hat * A1'
K_G_2 = A2 * K2hat * A2'

K_G = K_G_1 + K_G_2

% external reactions
Q = [Q1; Q2; Q3];

q = K_G\Q

D1 = A1' * q
D2 = A2' * q

d1 = lambda1 * D1
d2 = lambda2 * D2

f1 = K1 * d1
f2 = K2 * d2

F1 = lambda1' * f1
F2 = lambda2' * f2

figure(1)
plotDeflectedShape(0, 10, 0, 0, d1, magFactor, n, L, -90)
hold on
plotDeflectedShape(0, 0, 10, 0, d2, magFactor, n, L, 0)

figure(2)
plotDeflectedShapeModified(0, 10, d1, magFactor, n, L, -90)
hold on
plotDeflectedShapeModified(0, 0, d2, magFactor, n, L, 0)
