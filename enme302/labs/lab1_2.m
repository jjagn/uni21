% ENME302 lab 1 question 3
clear
close all
clc

% constants, material properties, values
E = 200 * 10^9; % elastic modulus, Pa
L1 = 1.1; % bar 1 length, m
L2 = 0.8; % bar 2 length, m
Q1 = 0;             % reactions, N
Q2 = -20 * 10^3;    % "" ""
A1 = 400; % bar 1 area, mm^2
A2 = 600; % bar 2 area, mm^2
alpha1 = 0; % angle of bar 1
alpha2 = 55; % angle of bar 2

% element 1
K1 = local_bar(E,A1,L1)
[K1hat, lambda1] = global_bar(K1, alpha1)

% element 2
K2 = local_bar(E,A2,L2)
[K2hat, lambda2] = global_bar(K2, alpha2)

% assembly matrices
A1 = [0 0 1 0; 0 0 0 1];
A2 = [1 0 0 0; 0 1 0 0];

% finding K_G
K_G_1 = A1 * K1hat * A1'
K_G_2 = A2 * K2hat * A2'

K_G_list = {K_G_1 K_G_2};

K_G = zeros(2);
for i = 1:length(K_G_list)
    K_G = K_G + K_G_list{i};
end

K_G

% external reactions
Q = [Q1; Q2];

q = K_G\Q

% finding nodal force vectors
F1 = K1hat * A1' * q
F2 = K2hat * A2' * q

% finding nodal displacement vectors
d1 = lambda1 * A1' * q
d2 = lambda2 * A2' * q

% finding strain
eps1 = (d1(1) - d1(2))/L1
eps2 = (d2(1) - d2(2))/(L2)

% finding nodal force vectors in local coordinates
f1 = K1 * d1
f2 = K2 * d2

% % plot system
% figure(1)
% 
% mag_factor = 1000; % magnification factor for displacements
% % original structure
% line([0, 10], [0, 10], 'Color', 'black')
% line([0, 10], [10, 10], 'Color', 'black')
% 
% % displaced structure
% title('FEA result')
% line([0, 10 + mag_factor * q(1)], [10, 10 + mag_factor * q(2)], 'Color', 'red', 'LineStyle', '-')
% line([0, 10 + mag_factor * q(1)], [0, 10 + mag_factor * q(2)], 'Color', 'red', 'LineStyle', '-')
% legend('original lmnt 1', 'original lmnt 2', 'disp. lmnt 1', 'disp. lmnt 2')