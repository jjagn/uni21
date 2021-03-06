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

L1 = 3;                  % frame length, m
L2 = L1;
L3 = L1;
L4 = L1;

A = 1*10^-3;            % frame area, m^2
I = 1*10^-5;            % frame area moment, m^4

alpha1 = 0;            % angle of frame 1
alpha2 = -75;            % angle of frame 2
% alpha3 = -20;
% alpha4 = -60;

magFactor = 10;
n = 50;
%=========================================================================%


%=========================================================================%
% CALCULATIONS
%=========================================================================%
% element 1
K1 = local_frame(E,I,A,L1);
[K1hat, lambda1] = global_frame(K1, alpha1);

% element 2
K2 = local_frame(E,I,A,L2);
[K2hat, lambda2] = global_frame(K2, alpha2);

% K3 = local_frame(E,I,A,L3);
% [K3hat, lambda3] = global_frame(K3, alpha3);
% 
% K4 = local_frame(E,I,A,L4);
% [K4hat, lambda4] = global_frame(K4, alpha4);
%=========================================================================%


%=========================================================================%
% ASSEMBLY MATRICES
A1 = eye(6)

A2 = [zeros(3, 6);
    eye(3), zeros(3)]

%=========================================================================%


%=========================================================================%
% FINDING K_G
K_G_1 = A1 * K1hat * A1';
K_G_2 = A2 * K2hat * A2';
% K_G_3 = A3 * K3hat * A3';
% K_G_4 = A4 * K4hat * A4';

K_G = K_G_1 + K_G_2;
%=========================================================================%

%=========================================================================%
% SOLVING SYSTEM FOR OVERALL DEFLECTIONS

f_eq_1 = addTransverseUDL(-2000, 3);

F_eq_1 = lambda1' * f_eq_1;

Q_UDL = A1 * F_eq_1;

Q_cable = 6e3 .*...
    [0;
     0;
     0;
     cosd(-35);
     sind(-35);
     0]

Qtotal = Q_UDL + Q_cable

q = K_G\Qtotal;

if printDeflections
    format long
    fprintf("q =\n")
    disp(q)
    format short
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying q, overall deflections\n")

end
%=========================================================================%


%=========================================================================%
% element nodal displacements, global coords

D1 = A1' * q;
D2 = A2' * q;
% D3 = A3' * q;
% D4 = A4' * q;

if printGlobalDisplacements
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("D, nodal displacements in global coords\n")
    fprintf("%s\n", spacerString)
    fprintf("D1 =\n")
    disp(D1)
    fprintf("D2 =\n")
    disp(D2)
%     fprintf("D3 =\n")
%     disp(D3)
%     fprintf("D4 =\n")
%     disp(D4)

else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying D, nodal displacements in global coords\n")

end
%=========================================================================%


%=========================================================================%
% element nodal displacements, local coords

d1 = lambda1 * D1;
d2 = lambda2 * D2;
% d3 = lambda3 * D3;
% d4 = lambda4 * D4;

if printElementDisplacements
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("d, nodal displacements in element coords\n")
    fprintf("%s\n", spacerString)
    
    fprintf("d1 =\n")
    disp(d1)
    fprintf("d2 =\n")
    disp(d2)
%     fprintf("d3 =\n")
%     disp(d3)
%     fprintf("d4 =\n")
%     disp(d4)
    

else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying d, nodal displacements in element coords\n")

end
%=========================================================================%


%=========================================================================%
% Element forcing vectors

f1 = K1 * d1;
f2 = K2 * d2;
% f3 = K3 * d3;
% f4 = K4 * d4;

if printElementForcingVectors
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("f, forcing vectors in element coords\n")
    fprintf("%s\n", spacerString)
    
    fprintf("f1 =\n")
    disp(f1)
    fprintf("f2 =\n")
    disp(f2)
%     fprintf("f3 =\n")
%     disp(f3)
%     fprintf("f4 =\n")
%     disp(f4)


else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying f, element forcing vectors in element coords\n")

end
%=========================================================================%


%=========================================================================%
% Global forcing vectors

F1 = lambda1' * f1;
F2 = lambda2' * f2;
% F3 = lambda3' * f3;
% F4 = lambda4' * f4;

if printGlobalForcingVectors
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("F, forcing vectors in global coords\n")
    fprintf("%s\n", spacerString)
    
    fprintf("F1 =\n")
    disp(F1)
    fprintf("F2 =\n")
    disp(F2)
%     fprintf("F3 =\n")
%     disp(F3)
%     fprintf("F4 =\n")
%     disp(F4)


else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying F, element forcing vectors in global coords\n")

end
%=========================================================================%


%=========================================================================%
% axial stress
if printStress
    fprintf("\n")
    fprintf("%s\n", spacerString)
    m1AxialLoadkN = f2(1) / 1000;
    m1AxialStress = f2(1) / A;
    m1AxialStressMPa = m1AxialStress *1e-6;
    if m1AxialLoadkN > 0
        fprintf("member 2 in compression, %.4fkN\n", m1AxialLoadkN)
        fprintf("axial stress = %.4fMPa\n", m1AxialStressMPa)
    elseif m1AxialLoadkN < 0 
        fprintf("member 2 in tension, %.4fkN\n", m1AxialLoadkN)
        fprintf("axial stress = %.4fMPa\n", m1AxialStressMPa)
    else 
        fprintf("member 2 under no axial load\n")
    end
    fprintf("\n")
    fprintf("%s\n", spacerString)
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying stresses\n")

end
%=========================================================================%


%=========================================================================%
% axial strain
if printStrain
    eps1 = (d1(1) - d1(4))/L1
    eps2 = (d2(1) - d2(4))/L2
%     eps3 = (d3(1) - d3(4))/L3
%     eps4 = (d4(1) - d4(4))/L4
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

    e2originX = e1originX + L1*cosd(alpha1);
    e2originY = e1originY + L1*sind(alpha1);

%     e3originX = e2originX + L2*cosd(alpha2);
%     e3originY = e2originY + L2*sind(alpha2);
% 
%     e4originX = e3originX + L3*cosd(alpha3);
%     e4originY = e3originY + L3*sind(alpha3);

    % plotting sketch
    plotDeflectedShapeModified(e1originX, e1originY, d1, magFactor, n, L1, alpha1)
    hold on
    plotDeflectedShapeModified(e2originX, e2originY, d2, magFactor, n, L2, alpha2)
%     plotDeflectedShapeModified(e3originX, e3originY, d3, magFactor, n, L3, alpha3)
%     plotDeflectedShapeModified(e4originX, e4originY, d4, magFactor, n, L4, alpha4)

    axis equal
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not plotting sketch\n")
    fprintf("%s\n", spacerString)
end
% ==========================================================
