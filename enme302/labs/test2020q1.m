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
printElementDisplacements =     0;     % d
printGlobalDisplacements =      0;     % D

printStrain =                   0;
printStress =                   0;
plotSketch =                    1;

spacerString = "==========================================================";
%=========================================================================%


%=========================================================================%
% INPUT VARIABLES
E = 10 * 10^9;         % elastic modulus, Pa

L1 = 0.5;    % bar 1 length, m
L2 = 2;                % bar 2 length, m
L3 = 2;                % bar 3 length, m
% L4 = 1554.92 / 1000;    % bar 4 length, m
% L5 = 1041.36 / 1000;    % bar 5 length, m
% L6 = 777.460 / 1000;    % bar 6 length, m
% L7 = 2;                 % bar 7 length, m
% L8 = 2;                 % bar 8 length, m

Ls = L1;
Ls(2) = L2;
Ls(3) = L3;
% Ls(4) = L4;
% Ls(5) = L5;
% Ls(6) = L6;
% Ls(7) = L7;
% Ls(8) = L8;

% D = 0.1;                        % major diameter of hollow circular bar
% d = 0.08;                       % minor diameter of hollow circular bar
% I = pi*(D^4-d^4)/64;            % bar area moment, m^4
% A = (pi*D^2)/4-(pi*d^2)/4;      % bar area, m^2

A = 4*10^-4;            % frame area, m^2
I = 1*10^-6;            % frame area moment, m^4

alpha1 = 90;         % angle of frame 1
alpha2 = 45;        % angle of frame 2
alpha3 = -75;        % angle of frame 3
% alpha4 = 30.97;         % 
% alpha5 = -50.18;        %
% alpha6 = 30.99;         %
% alpha7 = 0;             %
% alpha8 = 0;             %

alphas = alpha1;
alphas(2) = alpha2;
alphas(3) = alpha3;
% alphas(4) = alpha4;
% alphas(5) = alpha5;
% alphas(6) = alpha6;
% alphas(7) = alpha7;
% alphas(8) = alpha8;

n_elements = 3;
magFactor = 1;
n = 50;

%=========================================================================%


%=========================================================================%
% K_G CALCULATIONS
%=========================================================================%

for i = 1:n_elements
Ks(:,:,i) = local_frame(E,I,A,Ls(i));
[Khats(:,:,i), lambdas(:,:,i)] = global_frame(Ks(:,:,i), alphas(i));
end

%=========================================================================%


%=========================================================================%
% ASSEMBLY MATRICES
A1 = [zeros(3) eye(3); zeros(3,6)]

A2 = [eye(6)]

A3 = [zeros(3,6); eye(3) zeros(3)]
  
AssemblyMatrices = A1;
AssemblyMatrices(:,:,2) = A2;
AssemblyMatrices(:,:,3) = A3;
% AssemblyMatrices(:,:,4) = A4;
% AssemblyMatrices(:,:,5) = A5;
% AssemblyMatrices(:,:,6) = A6;
% AssemblyMatrices(:,:,7) = A7;
% AssemblyMatrices(:,:,8) = A8;
%=========================================================================%


%=========================================================================%
% FINDING K_G

K_G = 0;
for i = 1:n_elements
K_Gs(:,:,i) = AssemblyMatrices(:,:,i) * Khats(:,:,i) * AssemblyMatrices(:,:,i)';
K_G = K_G + K_Gs(:,:,i);
end

%=========================================================================%

%=========================================================================%
% SOLVING SYSTEM FOR OVERALL DEFLECTIONS
alpha_load = -80;
load_intensity = 4000;
transverse_component = load_intensity *sind(alpha_load)
axial_component = load_intensity * cosd(alpha_load)

f_eq_1 = addTransversePointLoad(transverse_component, L2, 1.5)

f_eq_2 = addAxialPointLoad(axial_component, L2, 1.5)

F_eq_1 = lambdas(:,:,2)' * f_eq_1;
F_eq_2 = lambdas(:,:,2)' * f_eq_2;

Q = A2 * F_eq_1 + A2 * F_eq_2

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
    stresses = zeros(1,8);
    fprintf("\n")
    fprintf("%s\n", spacerString)
    sum = 0;
    for i = 1:n_elements
        f = fs(:,:,i);
        AxialLoadkN = f(1)/ 1000;
        AxialStress = f(1) / A; % DID NOT HAVE TO VECTORISE A AS ALL HAVE SAME AREA
        AxialStressMPa = AxialStress *1e-6;
        if AxialLoadkN > 0
            fprintf("member %d in compression, %.6fkN\n", i, AxialLoadkN)
            fprintf("axial stress = %.6fMPa\n", AxialStressMPa)
        elseif AxialLoadkN < 0 
            fprintf("member %d in tension, %.6fkN\n", i, AxialLoadkN)
            fprintf("axial stress = %.6fMPa\n", AxialStressMPa)
        else 
            fprintf("member %d under no axial load\n", i)
        end
        fprintf("\n")
        fprintf("%s\n", spacerString)
        stresses(i) = AxialStressMPa;
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
        L = Ls(i);
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
    e1originY = 0;

    e2originX = e1originX + L1*cosd(alpha1);
    e2originY = e1originY + L1*sind(alpha1);

    e3originX = e2originX + L2*cosd(alpha2);
    e3originY = e2originY + L2*sind(alpha2);
% 
%     e4originX = 0;
%     e4originY = 0;
%     
%     e5originX = e4originX + L4*cosd(alpha4);
%     e5originY = e4originY + L4*sind(alpha4);
%     
%     e6originX = e5originX + L5*cosd(alpha5);
%     e6originY = e5originY + L5*sind(alpha5);
%     
%     e7originX = 0;
%     e7originY = 0;
%     
%     e8originX = e7originX + L7*cosd(alpha7);
%     e8originY = e7originY + L7*sind(alpha7);

    Xorigins = e1originX;
    Xorigins(2) = e2originX;
    Xorigins(3) = e3originX;
%     Xorigins(4) = e4originX;
%     Xorigins(5) = e5originX;
%     Xorigins(6) = e6originX;
%     Xorigins(7) = e7originX;
%     Xorigins(8) = e8originX;

    Yorigins = e1originY;
    Yorigins(2) = e2originY;
    Yorigins(3) = e3originY;
%     Yorigins(4) = e4originY;
%     Yorigins(5) = e5originY;
%     Yorigins(6) = e6originY;
%     Yorigins(7) = e7originY;
%     Yorigins(8) = e8originY;

    % plotting sketch
    for i = 1:n_elements
%         plotDeflectedBar(Xorigins(i), Yorigins(i), Ds(:,:,i), magFactor, n, Ls(i), alphas(i))
        plotDeflectedShapeModified(Xorigins(i), Yorigins(i), ds(:,:,i), magFactor, n, Ls(i), alphas(i))
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