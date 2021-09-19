clear
clc
close all
format short

% -------------------------- GLOBAL CONSTANTS ----------------------------
E = 200e9;         % Elastic modulus of the steel members in Pa
diam1 = 100/1000;   % Outer diameter of the member in m
diam2 = 80/1000;   % Inner diameter of the member in m
A = pi*(diam1^2 - diam2^2)/4;  % Cross-sectional area of the member in m^2
I = pi*(diam1^4 - diam2^4)/64; % Second area moment of members in m^4
n = 8;   % Number of Elements in structure

% -------------------------- ELEMENT CONSTANTS ---------------------------
% Element Lengths in m
L1 = sqrt((4/3)^2 + 0.4^2);
L2 = L1;
L3 = L1;
L4 = sqrt((4/3)^2 + 0.8^2);
L5 = sqrt((4/6)^2 + 0.8^2);
L6 = L4/2;
L7 = 2;
L8 = L7;

% Element Angles of Rotation in degrees
alpha1 = -atand(0.4/(4/3));
alpha2 = alpha1;
alpha3 = alpha1;
alpha4 = atand(0.8/(4/3));
alpha5 = -atand(0.8/(4/6));
alpha6 = alpha4;
alpha7 = 0;
alpha8 = 0;

% Element External Forcing Terms in N
Q = zeros(12, 1);
Q8 = -25000;  % The only externally applied load acts down at element 8
Q(8) = Q8;

% ---------------- STIFFNESS MATRICES IN ELEMENT COORDS ------------------
K1 = local_frame(E, I, A, L1);
K2 = local_frame(E, I, A, L2);
K3 = local_frame(E, I, A, L3);
K4 = local_frame(E, I, A, L4);
K5 = local_frame(E, I, A, L5);
K6 = local_frame(E, I, A, L6);
K7 = local_frame(E, I, A, L7);
K8 = local_frame(E, I, A, L8);


% ------------- ELEMENT STIFFNESS MATRICES IN GLOBAL COORDS --------------
[K1hat, Lambda1] = global_frame(K1, alpha1);
[K2hat, Lambda2] = global_frame(K2, alpha2);
[K3hat, Lambda3] = global_frame(K3, alpha3);
[K4hat, Lambda4] = global_frame(K4, alpha4);
[K5hat, Lambda5] = global_frame(K5, alpha5);
[K6hat, Lambda6] = global_frame(K6, alpha6);
[K7hat, Lambda7] = global_frame(K7, alpha7);
[K8hat, Lambda8] = global_frame(K8, alpha8);

% ---------------------------- ASSEMBLY MATRICES -------------------------
zeroM = zeros(12,6); % Zero matrix for 12 global DOF's % 6 Element DOF's

A1 = zeroM; % Initialising Element 1 assembly matrix
A1(1:3, 4:6) = eye(3);

A2 = zeroM; % Initialising Element 2 assembly matrix
A2(1:6, 1:6) = eye(6);

A3 = zeroM; % Initialising Element 3 assembly matrix
A3(4:9, 1:6) = eye(6);

A4 = zeroM; % Initialising Element 4 assembly matrix
A4(1:3, 4:6) = eye(3);

A5 = zeroM; % Initialising Element 5 assembly matrix
A5(1:3, 1:3) = eye(3);
A5(10:12, 4:6) = eye(3);

A6 = zeroM; % Initialising Element 6 assembly matrix
A6(4:6, 4:6) = eye(3);
A6(10:12, 1:3) = eye(3);

A7 = zeroM; % Initialising Element 7 assembly matrix
A7(10:12, 4:6) = eye(3);

A8 = zeroM; % Initialising Element 8 assembly matrix
A8(7:9, 4:6) = eye(3);
A8(10:12, 1:3) = eye(3);


% --------- ELEMENTS CONTRIBUTIONS TO THE GLOBAL STIFFNESS MATRIX --------
K1G = A1*K1hat*A1';
K2G = A2*K2hat*A2';
K3G = A3*K3hat*A3';
K4G = A4*K4hat*A4';
K5G = A5*K5hat*A5';
K6G = A6*K6hat*A6';
K7G = A7*K7hat*A7';
K8G = A8*K8hat*A8';

% --------------------- TOTAL GLOBAL STIFFNESS MATRIX --------------------
KG = K1G + K2G + K3G + K4G + K5G + K6G + K7G + K8G;

% --------------------- OVERALL STRUCTURAL DEFLECTIONS -------------------
q = KG\Q;       % ---- SOLVED -----


% --------------------------- POST PROCESSING ----
% Global deflection vectors
D1 = (A1')*q;  % Global deflection vector for element 1
D2 = (A2')*q;
D3 = (A3')*q;
D4 = (A4')*q;
D5 = (A5')*q;
D6 = (A6')*q;
D7 = (A7')*q;
D8 = (A8')*q;


% Local deflection vectors
d1 = Lambda1*D1;  % Local deflection vector for element 1
d2 = Lambda2*D2;
d3 = Lambda3*D3;
d4 = Lambda4*D4;
d5 = Lambda5*D5;
d6 = Lambda6*D6;
d7 = Lambda7*D7;
d8 = Lambda8*D8;


% Element nodal forcing vectors in element coords
f1 = K1*d1;  % Element nodal forcing vector for element 1 - element coords
f2 = K2*d2;
f3 = K3*d3;
f4 = K4*d4;
f5 = K5*d5;
f6 = K6*d6;
f7 = K7*d7;
f8 = K8*d8;

f = [f1, f2, f3, f4, f5, f6, f7, f8];

% Element nodal forcing vectors in global coords
F1 = K1hat*D1; % Element nodal forcing vector for element 1 - global coords
F2 = K2hat*D2;
F3 = K3hat*D3;
F4 = K4hat*D4;
F5 = K5hat*D5;
F6 = K6hat*D6;
F7 = K7hat*D7;
F8 = K8hat*D8;


% ------------------ DEFLECTIONS AT RIGHT HAND TIP -----------------------
tipDeflectionX = q(7)*1000  % Deflection of tip in x direction in mm
tipDeflectionY = q(8)*1000  % Deflection of tip in x direction in mm
tipDeflectionRot = q(9)*pi*1000/180  % Rotational deflection in mrad


% ------------------ REACTION FORCES AT SUPPORT A ------------------------
FAx = F1(1)/1000   % XG-component of reaction force at support A
FAy = F1(2)/1000   % YG-component of reaction force at support A
MA = F1(3)/1000    % Reaction moment at support A


% ------------------ REACTION FORCES AT SUPPORT B ------------------------
FBx = (F4(1) + F7(1))/1000 % XG-component of reaction force at support B
FBy = (F4(2) + F7(2))/1000 % YG-component of reaction force at support B
MB = (F4(3) + F7(3))/1000  % Reaction moment at support B


% ---------------------------- AXIAL FORCES ------------------------------
axialForceAbs1 = abs(f1(1)); % Abs value of the axial force in element 1
axialForceAbs2 = abs(f2(1)); % Abs value of the axial force in element 2
axialForceAbs3 = abs(f3(1)); % Abs value of the axial force in element 3
axialForceAbs4 = abs(f4(1)); % Abs value of the axial force in element 4
axialForceAbs5 = abs(f5(1)); % Abs value of the axial force in element 5
axialForceAbs6 = abs(f6(1)); % Abs value of the axial force in element 6
axialForceAbs7 = abs(f7(1)); % Abs value of the axial force in element 7
axialForceAbs8 = abs(f8(1)); % Abs value of the axial force in element 8

% --------------------- AXIAL STRESSES AND STRAINS -----------------------
axialStrain1 = (d1(4) - d1(1))/L1;  % Axial Strain of element 1
axialStress1 = E*axialStrain1;     % Axial Stress of element 1 in Pa
axialStressMPa1 = axialStress1/(1e6); % Axial Stress of element 1 in MPa

axialStrain2 = (d2(4) - d2(1))/L2;  % Axial Strain of element 2
axialStress2 = E*axialStrain2;      % Axial Stress of element 2 in Pa
axialStressMPa2 = axialStress2/(1e6); % Axial Stress of element 2 in MPa

axialStrain3 = (d3(4) - d3(1))/L3;  % Axial Strain of element 3
axialStress3 = E*axialStrain3;      % Axial Stress of element 3 in Pa
axialStressMPa3 = axialStress3/(1e6); % Axial Stress of element 3 in MPa

axialStrain4 = (d4(4) - d4(1))/L4;  % Axial Strain of element 4
axialStress4 = E*axialStrain4;      % Axial Stress of element 4 in Pa
axialStressMPa4 = axialStress4/(1e6); % Axial Stress of element 4 in MPa

axialStrain5 = (d5(4) - d5(1))/L5;  % Axial Strain of element 5
axialStress5 = E*axialStrain5;      % Axial Stress of element 5 in Pa
axialStressMPa5 = axialStress5/(1e6); % Axial Stress of element 5 in MPa

axialStrain6 = (d6(4) - d6(1))/L6;  % Axial Strain of element 6
axialStress6 = E*axialStrain6;      % Axial Stress of element 6 in Pa
axialStressMPa6 = axialStress6/(1e6); % Axial Stress of element 6 in MPa

axialStrain7 = (d7(4) - d7(1))/L7;  % Axial Strain of element 7
axialStress7 = E*axialStrain7;      % Axial Stress of element 7  in Pa
axialStressMPa7 = axialStress7/(1e6); % Axial Stress of element 7 in MPa

axialStrain8 = (d8(4) - d8(1))/L8;  % Axial Strain of element 8
axialStress8 = E*axialStrain8;      % Axial Stress of element 8 in Pa
axialStressMPa8 = axialStress8/(1e6); % Axial Stress of element 8 in MPa

% ------------------------ BENDING STRESSES ------------------------------
c = diam1/2;   % Distance to the outermost fibre for bending in m

bendingStress1 = abs(f1(3)*c/I)/(1e6); % Bending stress of element 1 in MPa
bendingStress2 = abs(f2(3)*c/I)/(1e6); % Bending stress of element 2 in MPa
bendingStress3 = abs(f3(3)*c/I)/(1e6); % Bending stress of element 3 in MPa
bendingStress4 = abs(f4(3)*c/I)/(1e6); % Bending stress of element 4 in MPa
bendingStress5 = abs(f5(3)*c/I)/(1e6); % Bending stress of element 5 in MPa
bendingStress6 = abs(f6(3)*c/I)/(1e6); % Bending stress of element 6 in MPa
bendingStress7 = abs(f7(3)*c/I)/(1e6); % Bending stress of element 7 in MPa
bendingStress8 = abs(f8(3)*c/I)/(1e6); % Bending stress of element 8 in MPa

% ------------------------- TOTAL STRESSES -------------------------------
totalStress1 = abs(axialStressMPa1) + bendingStress1;
totalStress2 = abs(axialStressMPa2) + bendingStress2;
totalStress3 = abs(axialStressMPa3) + bendingStress3;
totalStress4 = abs(axialStressMPa4) + bendingStress4;
totalStress5 = abs(axialStressMPa5) + bendingStress5;
totalStress6 = abs(axialStressMPa6) + bendingStress6;
totalStress7 = abs(axialStressMPa7) + bendingStress7;
totalStress8 = abs(axialStressMPa8) + bendingStress8;

Max_Stress = max([totalStress1, totalStress2, totalStress3, totalStress4, totalStress5,...
    totalStress5, totalStress7, totalStress8])

delta = [tipDeflectionX, tipDeflectionY, tipDeflectionRot, FAx, FAy, MA, FBx, FBy, MB, Max_Stress]