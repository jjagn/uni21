M = 1.5; % mass of cart, kg
km = 0.017; % motor back emf constant, V/rad/s
kg = 3.7; % gearing ratio
R = 1.5; % resistance of motor armature, ohms
r = 0.018; % radius of pinion, m
D = 7; % damping present in physical system, e.g friction

% constants so I don't have to type all that junk out
B = (km * kg) / (m * R * r);
C = D / M + (km^2 * kg^2) / (m * R * r^2);

% controller gains
Kd = 1;
Kp = 100;

% fraction declaration
num = [B*Kd B*Kp];
den = [1 C+B*Kd B*Kp];

step(num, den);