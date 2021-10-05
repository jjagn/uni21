clear; close all; clc

% VARIABLES
L = 2;
H = 3;

f = @(x, y) x .* y .* (L-x) .* (H-y);
g = 0;
s = 0;
c = 6;
n_space = 30; % spatial grid points
tn = 300; % number of time points, e.g frames
desired_fps = 240
t_f = tn/desired_fps; % end time of sim

% ANALYTICAL SOLUTION
a_lim = 5;      % number of fourier series terms
b_lim = 5;      % "" ""

u_a = @(x, y, t) 0;


for b = 1:1:b_lim
    for a = 1:1:a_lim
        u_a = @(x, y, t) u_a(x, y, t) + (576/pi^6) .* (((1 + (-1)^(a+1)) * (1+(-1)^(b+1))) / (a^3*b^3)).*sin(a.*pi./2.*x).*sin(b.*pi./3.*y).*cos(pi.*sqrt(9.*a.^2+4.*b.^2).*t);
    end
end

x_points = linspace(0,L,n_space);
y_points = linspace(0,H,n_space);
t_points = linspace(0,t_f,tn);

[X,Y] = meshgrid(x_points, y_points);

for i = 1:tn
    results(:,:,i) = u_a(X, Y, t_points(i));
end

figure()
filename = 'analytical_soln_better_fps_2.gif';

for n = 1:tn
    surf(results(:,:,n))
    axis([0 n_space 0 n_space -2.5 2.5]);
    % pause(1/60)
    if n == 1
        gif(filename)
    else 
        gif('DelayTime', 1/60)
    end
    clc
    fprintf("iteration %d of %d\n", n, tn)
    fprintf("progress: %.2f%% \n", n/tn*100)
end
