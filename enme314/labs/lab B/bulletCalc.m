% Calculating a bullet's trajectory given initial velocity, drag
% coefficient, and starting height, for ENME31421S1 Lab B.
%
% 09.05.21, Jackson Crawford


clear
close all
clc

% assuming constant air density, bullet is fired ~at sea level, ground
% stays flat and level for the bullet's full flight. Also assumes that the
% drag coefficient stays constant, even though the reynolds number will
% change throughout the bullet's flight.


% bullet parameters
Cd = 0.49; % drag coefficient for bullet, obtained from sphere results and similarity 0.49
d = 5 * 10^-3; % bullet diameter, m
A = (pi * d^2)/4; % cross sectional area of bullet, m^2
m = 0.5*10^-3; % 0.5 grams in kg
mu = 1.519 * 10^-5;

v0 = 300; % initial velocity of bullet, ms^-1
x(1) = 0; % initial position in x, m
y(1) = 1.8; % initial position in y, m (bullet fired 1.8m from ground)
theta = 0; % angle at which bullet was initially fired at, degrees
vx = v0*cosd(theta); % initial velocity in x, ms^-1
vy = v0*sind(theta); % initial velocity in y, ms^-1

% environment parameters
rho = 1.207; % density of  air at sea level, kgm^-3
g = -9.81; % acceleration due to gravity, ms^-2

% program parameters

dt = 0.001; % time step, s
t(1) = 0; % initial time

% iterative solution

i = 1;
while min(y) >= 0
    
    Re = (rho * d * sqrt(vx^2 + vy^2))/mu;
    
    Cd = 24/Re * (1 + 0.15 * Re^0.687) + 0.42/(1 + 4.25 * 10^4 * Re^(-1.16));
    
    v_dot_x = - ((rho * Cd * A/2) * vx^2)/m;
    v_dot_y = g - ((rho * Cd * A/2) * vy^2)/m ;
    
    vx = vx + v_dot_x * dt;
    vy = vy + v_dot_y * dt;
    
    x(i+1) = x(i) + vx * dt + 1/2 * v_dot_x * dt^2;
    y(i+1) = y(i) + vy * dt + 1/2 * v_dot_y * dt^2;
    t(i+1) = t(i) + dt;
    
    i = i + 1;
end

% plotting
figure('DefaultAxesFontSize', 14)
plot(t, y)
ylabel('Bullet height (m)', 'Fontname', 'SF Pro Text')
xlabel('Time elapsed (s)', 'Fontname', 'SF Pro Text')
ylim([0, 1.85])

figure('DefaultAxesFontSize', 14)
plot(x, y)
ylabel('Bullet altitude (m)', 'Fontname', 'SF Pro Text')
xlabel('Traveled distance (m)', 'Fontname', 'SF Pro Text')
ylim([0, 1.85])
x(i)
y(i)
t(i)