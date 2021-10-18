clear; close all; clc



% VARIABLES
L = 2;
H = 3;

u_e = @(x, y, t) x.*(L-x).*y.*(H-y).*(1+1/2.*t);

g = @(x, y, t) 0.5 .* u_e;
s = @(x, y, t) 2.*c^2 .* (1 + 0.5.*t) .* (y .* (H-y) + x .* (L-x)); 
c = 6;



lambda = 0.7;
delta = 0.1;
grid_size = 1; % physical size of grid
t_f = 10; % end time of sim

n = round(grid_size/delta); % calculating spatial grid points to maintain constant courant number

dt = lambda*delta/c;
t_n = t_f/dt;


x = linspace(0,L,L/delta+1);
y = linspace(0,H,H/delta+1);
t = linspace(0,t_f,t_n);

[X,Y] = meshgrid(x, y);

for i = 1:t_n
    results_analytical(:,:,i) = u_e(X, Y, t(i));
    surf(results_analytical(:,:,i));
    pause(1/60)
end

