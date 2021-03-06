% jackson crawford
% enme303 assignment 1

clear
clc
close all

% PD control algorithm using syms and laplace transforms to keep track of
% voltage

% declaring symbolic variables
syms Kp Kd e x;

% function for voltage and laplace transform
v = Kp*e + Kd*diff(e);

V = laplace(v);

% cart dynamics

M = 1.5; % mass of cart, kg
Km = 0.017; % motor back emf constant, V/rad/s
Kg = 3.7; % gearing ratio
R = 1.5; % resistance of motor armature, ohms
r = 0.018; % radius of pinion, m
D = 7; % damping present in physical system, e.g friction

% constants so I don't have to type all that junk out
B = (Km * Kg) / (M * R * r);
C = D / M + (Km^2 * Kg^2) / (M * R * r^2);

eqn1 = M*diff(x, 2) == B*v - C*diff(x)