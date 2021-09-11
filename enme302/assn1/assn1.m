%=========================================================================%
% ENME302 assignment 1 - Jackson Crawford
clear
close all
clc
%=========================================================================%


%=========================================================================%
% settings
printDeflections =              1;     % q
printElementForcingVectors =    0;     % f
printGlobalForcingVectors =     0;     % F
printElementDisplacements =     1;     % d
printGlobalDisplacements =      0;     % D

printStrain =                   1;
printStress =                   1;
plotSketch =                    1;

spacerString = "==========================================================";
%=========================================================================%


%=========================================================================%
% INPUT VARIABLES
E = 200 * 10^9;         % elastic modulus, Pa

L1 = 1392.04 / 1000;    % bar 1 length, m
L2 = L1;                % bar 2 length, m
L3 = L1;                % bar 3 length, m
L4 = 1554.92 / 1000;    % bar 4 length, m
L5 = 1041.36 / 1000;    % bar 5 length, m
L6 = 777.460 / 1000;    % bar 6 length, m
L7 = 2;                 % bar 7 length, m
L8 = 2;                 % bar 8 length, m

Ls = L1;
Ls(:,:,2) = L2;
Ls(:,:,3) = L3;
Ls(:,:,4) = L4;
Ls(:,:,5) = L5;
Ls(:,:,6) = L6;
Ls(:,:,7) = L7;
Ls(:,:,8) = L8;

D = 0.05;               % major diameter of hollow circular bar
d = 0.038;              % minor diameter of hollow circular bar
I = pi*(D^4-d^4)/64;    % bar area moment, m^4
A = (pi*D^2)/4-(pi*d^2)/4;           % bar area, m^2

alpha1 = -16.7;         % angle of frame 1
alpha2 = alpha1;        % angle of frame 2
alpha3 = alpha1;        % angle of frame 3
alpha4 = 30.97;         % 
alpha5 = -50.18;        %
alpha6 = 30.99;         %
alpha7 = 0;             %
alpha8 = 0;             %

alphas = alpha1;
alphas(:,:,2) = alpha2;
alphas(:,:,3) = alpha3;
alphas(:,:,4) = alpha4;
alphas(:,:,5) = alpha5;
alphas(:,:,6) = alpha6;
alphas(:,:,7) = alpha7;
alphas(:,:,8) = alpha8;


n_elements = 8;
magFactor = 30;
n = 50;

%=========================================================================%


%=========================================================================%
% CALCULATIONS
%=========================================================================%
% element 1
K1 = local_bar(E,A,L1);
[K1hat, lambda1] = global_bar(K1, alpha1);

% element 2
K2 = local_bar(E,A,L2);
[K2hat, lambda2] = global_bar(K2, alpha2);

% element 3
K3 = local_bar(E,A,L3);
[K3hat, lambda3] = global_bar(K3, alpha3);

% element 4
K4 = local_bar(E,A,L4);
[K4hat, lambda4] = global_bar(K4, alpha4);

% element 5
K5 = local_bar(E,A,L5);
[K5hat, lambda5] = global_bar(K5, alpha5);

% element 6
K6 = local_bar(E,A,L6);
[K6hat, lambda6] = global_bar(K6, alpha6);

% element 7
K7 = local_bar(E,A,L7);
[K7hat, lambda7] = global_bar(K7, alpha7);

% element 8
K8 = local_bar(E,A,L8);
[K8hat, lambda8] = global_bar(K8, alpha8);

Ks = K1;
Ks(:,:,2) = K2;
Ks(:,:,3) = K3;
Ks(:,:,4) = K4;
Ks(:,:,5) = K5;
Ks(:,:,6) = K6;
Ks(:,:,7) = K7;
Ks(:,:,8) = K8;

lambdas = lambda1;
lambdas(:,:,2) = lambda2;
lambdas(:,:,3) = lambda3;
lambdas(:,:,4) = lambda4;
lambdas(:,:,5) = lambda5;
lambdas(:,:,6) = lambda6;
lambdas(:,:,7) = lambda7;
lambdas(:,:,8) = lambda8;
%=========================================================================%


%=========================================================================%
% ASSEMBLY MATRICES
A1 = [zeros(2) eye(2);
      zeros(6, 4)];

A2 = [eye(4);
      zeros(4)];
  
A3 = [zeros(2, 4);
      eye(4);
      zeros(2, 4)];
  
A4 = A1;

A5 = [eye(2) zeros(2);
      zeros(4, 4);
      zeros(2) eye(2)];
  
A6 = [zeros(2, 4);
      zeros(2) eye(2);
      zeros(2, 4);
      eye(2) zeros(2)];
  
A7 = [zeros(6, 4);
      zeros(2) eye(2)];
  
A8 = [zeros(4);
      zeros(2) eye(2);
      eye(2) zeros(2)];
  
AssemblyMatrices = A1;
AssemblyMatrices(:,:,2) = A2;
AssemblyMatrices(:,:,3) = A3;
AssemblyMatrices(:,:,4) = A4;
AssemblyMatrices(:,:,5) = A5;
AssemblyMatrices(:,:,6) = A6;
AssemblyMatrices(:,:,7) = A7;
AssemblyMatrices(:,:,8) = A8;
%=========================================================================%


%=========================================================================%
% FINDING K_G
K_G_1 = A1 * K1hat * A1';
K_G_2 = A2 * K2hat * A2';
K_G_3 = A3 * K3hat * A3';
K_G_4 = A4 * K4hat * A4';
K_G_5 = A5 * K5hat * A5';
K_G_6 = A6 * K6hat * A6';
K_G_7 = A7 * K7hat * A7';
K_G_8 = A8 * K8hat * A8';

K_G = K_G_1 + K_G_2 + K_G_3 + K_G_4 + K_G_5 + K_G_6 + K_G_7 + K_G_8;
%=========================================================================%

%=========================================================================%
% SOLVING SYSTEM FOR OVERALL DEFLECTIONS

Q = [0;...
     0;...
     0;...
     0;...
     0;...
     -25000;...
     0;...
     0];

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

end
%=========================================================================%


%=========================================================================%
% element nodal displacements, global coords
for i = 1:n_elements
Ds(:,:,i) = AssemblyMatrices(:,:,i)' * q;
end

if printGlobalDisplacements
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("D, nodal displacements in global coords\n")
    fprintf("%s\n", spacerString)
    
    for i = 1:n_elements
        fprintf("D%d =\n", i)
        disp(Ds(:,:,i))
    end

else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying D, nodal displacements in global coords\n")

end
%=========================================================================%


%=========================================================================%
% element nodal displacements, local coords

% d1 = lambda1 * D1;
% d2 = lambda2 * D2;
% d3 = lambda3 * D3;
% d4 = lambda4 * D4;
% d5 = lambda5 * D5;
% d6 = lambda6 * D6;
% d7 = lambda7 * D7;
% d8 = lambda8 * D8;

for i = 1:n_elements
ds(:,:,i) = lambdas(:,:,i) * Ds(:,:,i);
end

if printElementDisplacements
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("d, nodal displacements in element coords\n")
    fprintf("%s\n", spacerString)
    for i = 1:n_elements
        fprintf("d%d =\n", i)
        disp(ds(:,:,i))
    end
    
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying d, nodal displacements in element coords\n")

end
%=========================================================================%


%=========================================================================%
% Element forcing vectors
% 
% f1 = K1 * d1;
% f2 = K2 * d2;
% f3 = K3 * d3;
% f4 = K4 * d4;
% f5 = K5 * d5;
% f6 = K6 * d6;
% f7 = K7 * d7;
% f8 = K8 * d8;

for i = 1:n_elements
fs(:,:,i) = Ks(:,:,i) * ds(:,:,i);
end

if printElementForcingVectors
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("f, forcing vectors in element coords\n")
    fprintf("%s\n", spacerString)
    
    for i = 1:n_elements
        fprintf("f%d =\n", i)
        disp(fs(:,:,i))
    end


else
    
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying f, element forcing vectors in element coords\n")

end
%=========================================================================%


%=========================================================================%
% Global forcing vectors
for i = 1:n_elements
Fs(:,:,i) = lambdas(:,:,i)' * fs(:,:,i);
end

if printGlobalForcingVectors
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("F, forcing vectors in global coords\n")
    fprintf("%s\n", spacerString)
    
    for i = 1:n_elements
        fprintf("F%d =\n", i)
        disp(Fs(:,:,i))
    end


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
    for i = 1:n_elements
        f = fs(:,:,i);
        AxialLoadkN = f(1)/ 1000;
        AxialStress = f(1) / A; % DID NOT HAVE TO VECTORISE A AS ALL HAVE SAME AREA
        AxialStressMPa = AxialStress *1e-6;
        if AxialLoadkN > 0
            fprintf("member %d in compression, %.4fkN\n", i, AxialLoadkN)
            fprintf("axial stress = %.4fMPa\n", AxialStressMPa)
        elseif AxialLoadkN < 0 
            fprintf("member %d in tension, %.4fkN\n", i, AxialLoadkN)
            fprintf("axial stress = %.4fMPa\n", AxialStressMPa)
        else 
            fprintf("member %d under no axial load\n", i)
        end
        fprintf("\n")
        fprintf("%s\n", spacerString)
    end
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not displaying stresses\n")

end
%=========================================================================%


%=========================================================================%
% axial strain
if printStrain
    
    epss = zeros(1, 1, n_elements);
    for i = 1:n_elements
        d = ds(:,:,i);
        L = Ls(:,:,i);
        epss(:,:,i) = d(1) - d(2)/L;
    end
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
    e1originY = 1.2;

    e2originX = e1originX + L1*cosd(alpha1);
    e2originY = e1originY + L1*sind(alpha1);

    e3originX = e2originX + L2*cosd(alpha2);
    e3originY = e2originY + L2*sind(alpha2);

    e4originX = 0;
    e4originY = 0;
    
    e5originX = e4originX + L4*cosd(alpha4);
    e5originY = e4originY + L4*sind(alpha4);
    
    e6originX = e5originX + L5*cosd(alpha5);
    e6originY = e5originY + L5*sind(alpha5);
    
    e7originX = 0;
    e7originY = 0;
    
    e8originX = e7originX + L7*cosd(alpha7);
    e8originY = e7originY + L7*sind(alpha7);

    Xorigins = e1originX;
    Xorigins(:,:,2) = e2originX;
    Xorigins(:,:,3) = e3originX;
    Xorigins(:,:,4) = e4originX;
    Xorigins(:,:,5) = e5originX;
    Xorigins(:,:,6) = e6originX;
    Xorigins(:,:,7) = e7originX;
    Xorigins(:,:,8) = e8originX;

    Yorigins = e1originY;
    Yorigins(:,:,2) = e2originY;
    Yorigins(:,:,3) = e3originY;
    Yorigins(:,:,4) = e4originY;
    Yorigins(:,:,5) = e5originY;
    Yorigins(:,:,6) = e6originY;
    Yorigins(:,:,7) = e7originY;
    Yorigins(:,:,8) = e8originY;

    % plotting sketch
    for i = 1:n_elements
        plotDeflectedBar(Xorigins(:,:,i), Yorigins(:,:,i), Ds(:,:,i), magFactor, n, Ls(:,:,i), alphas(:,:,i))
        hold on
    end
    
    axis equal
else
    fprintf("\n")
    fprintf("%s\n", spacerString)
    fprintf("not plotting sketch\n")
    fprintf("%s\n", spacerString)
end
% ==========================================================

