%=========================================================================%
% ENME302
clear
close all
clc
%=========================================================================%


%=========================================================================%
% settings
printDeflections =              1;     % q
printElementForcingVectors =    1;     % f
printGlobalForcingVectors =     1;     % F
printElementDisplacements =     1;     % d
printGlobalDisplacements =      1;     % D

printStrain =                   1;
printStress =                   1;
plotSketch =                    1;

spacerString = "==========================================================";
%=========================================================================%


%=========================================================================%
% INPUT VARIABLES
E = 200 * 10^9;         % elastic modulus, Pa
L = 5;                  % frame length, m
A = 6*10^-4;            % frame area, m^2
I = 8*10^-6;            % frame area moment, m^4

alpha1 = 60;            % angle of frame 1
alpha2 = 20;            % angle of frame 2
alpha3 = -20;
alpha4 = -60;

magFactor = 10;
n = 50;
%=========================================================================%


%=========================================================================%
% CALCULATIONS
%=========================================================================%
% element 1
K1 = local_frame(E,I,A,L);
[K1hat, lambda1] = global_frame(K1, alpha1);

% element 2
K2 = local_frame(E,I,A,L);
[K2hat, lambda2] = global_frame(K2, alpha2);

K3 = local_frame(E,I,A,L);
[K3hat, lambda3] = global_frame(K3, alpha3);

K4 = local_frame(E,I,A,L);
[K4hat, lambda4] = global_frame(K4, alpha4);
%=========================================================================%


%=========================================================================%
% ASSEMBLY MATRICES
A1 = [zeros(3) eye(3);
    zeros(6)];

A2 = [eye(6);
    zeros(3, 6)];

A3 = [zeros(3, 6);
    eye(6)];

A4 = [zeros(6);
    eye(3) zeros(3)];

Q = [0;
    0;
    0;
    0;
    -5e4;
    0;
    0;
    0;
    0];
%=========================================================================%


%=========================================================================%
% FINDING K_G
K_G_1 = A1 * K1hat * A1';
K_G_2 = A2 * K2hat * A2';
K_G_3 = A3 * K3hat * A3';
K_G_4 = A4 * K4hat * A4';

K_G = K_G_1 + K_G_2 + K_G_3 + K_G_4;
%=========================================================================%

%=========================================================================%
% SOLVING SYSTEM FOR OVERALL DEFLECTIONS
q = K_G\Q;

if printDeflections
    format long
    fprintf("q =\n")
    disp(q)
    format short
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying q, overall deflections\n")
    fprintf("%s\n", spacerString)
end
%=========================================================================%


%=========================================================================%
% element nodal displacements, global coords

D1 = A1' * q;
D2 = A2' * q;
D3 = A3' * q;
D4 = A4' * q;

if printGlobalDisplacements
    fprintf("D1 =\n")
    disp(D1)
    fprintf("D2 =\n")
    disp(D2)
    fprintf("D3 =\n")
    disp(D3)
    fprintf("D4 =\n")
    disp(D4)
else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying D, nodal displacements in global coords\n")
    fprintf("%s\n", spacerString)
end
%=========================================================================%


%=========================================================================%
% element nodal displacements, local coords

d1 = lambda1 * D1;
d2 = lambda2 * D2;
d3 = lambda3 * D3;
d4 = lambda4 * D4;

if printElementDisplacements
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("d, nodal displacements in element coords\n")
    fprintf("%s\n", spacerString)
    
    fprintf("d1 =\n")
    disp(d1)
    fprintf("d2 =\n")
    disp(d2)
    fprintf("d3 =\n")
    disp(d3)
    fprintf("d4 =\n")
    disp(d4)
    
    fprintf("%s\n", spacerString)
else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying d, nodal displacements in element coords\n")
    fprintf("%s\n", spacerString)
end
%=========================================================================%


%=========================================================================%
% Element forcing vectors

f1 = K1 * d1;
f2 = K2 * d2;
f3 = K3 * d3;
f4 = K4 * d4;

if printElementForcingVectors
    disp(f1)
    disp(f2)
    disp(f3)
    disp(f4)
else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying f, element forcing vectors in element coords\n")
    fprintf("%s\n", spacerString)
end
%=========================================================================%


%=========================================================================%
% Global forcing vectors

F1 = lambda1' * f1;
F2 = lambda2' * f2;
F3 = lambda3' * f3;
F4 = lambda4' * f4;

if printGlobalForcingVectors
    disp(F1)
    disp(F2)
    disp(F3)
    disp(F4)
else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying F, element forcing vectors in global coords\n")
    fprintf("%s\n", spacerString)
end
%=========================================================================%


%=========================================================================%
% axial stress
if printStress
    fprintf("\n")
    fprintf("%s\n", spacerString)
    m1AxialStresskN = f1(1) / 1000;
    if m1AxialStresskN > 0
        fprintf("member 1 in compression, %.4fkN", m1AxialStresskN)
    elseif m1AxialStresskN < 0 
        fprintf("member 1 in tension, %.4fkN", m1AxialStresskN)
    else 
        fprintf("member 1 under no axial load")
    end
    fprintf("\n")
    fprintf("%s\n", spacerString)
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying stresses\n")
    fprintf("%s\n", spacerString)
end
%=========================================================================%


%=========================================================================%
% axial strain
if printStrain
    eps1 = (d1(1) - d1(4))/L
    eps2 = (d2(1) - d2(4))/L
    eps3 = (d3(1) - d3(4))/L
    eps4 = (d4(1) - d4(4))/L
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying strains\n")
    fprintf("%s\n", spacerString)
end
% ==========================================================


% ==========================================================
if plotSketch
    % sketching calculations
    e1originX = 0;
    e1originY = 0;

    e2originX = e1originX + L*cosd(alpha1);
    e2originY = e1originY + L*sind(alpha1);

    e3originX = e2originX + L*cosd(alpha2);
    e3originY = e2originY + L*sind(alpha2);

    e4originX = e3originX + L*cosd(alpha3);
    e4originY = e3originY + L*sind(alpha3);

    % plotting sketch
    plotDeflectedShapeModified(e1originX, e1originY, d1, magFactor, n, L, alpha1)
    hold on
    plotDeflectedShapeModified(e2originX, e2originY, d2, magFactor, n, L, alpha2)
    plotDeflectedShapeModified(e3originX, e3originY, d3, magFactor, n, L, alpha3)
    plotDeflectedShapeModified(e4originX, e4originY, d4, magFactor, n, L, alpha4)

    axis equal
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not plotting sketch\n")
    fprintf("%s\n", spacerString)
end
% ==========================================================
