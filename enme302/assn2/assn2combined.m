clc; close all; clear

% VARIABLES
L = 2;
H = 3;

g = 0;
s = 0;
c = 6;
n = 100; % spatial grid points
grid_size = 1; % physical size of grid
t_n = 1000; % number of time points
t_f = 1; % end time of sim

animated_plot = 0;

delta = grid_size/n;
dt = t_f/t_n;

lambda = c*dt/delta
lambda2 = lambda*lambda

% ANALYTICAL SOLUTION
a_lim = 5;      % number of fourier series terms
b_lim = 5;      % "" ""

u_a = @(x, y, t) 0;

for b = 1:1:b_lim
    for a = 1:1:a_lim
        u_a = @(x, y, t) u_a(x, y, t) + (576/pi^6) .* (((1 + (-1)^(a+1)) * (1+(-1)^(b+1))) / (a^3*b^3)).*sin(a.*pi./2.*x).*sin(b.*pi./3.*y).*cos(pi.*sqrt(9.*a.^2+4.*b.^2).*t);
    end
end

x = linspace(0,L,L/delta+1);
y = linspace(0,H,H/delta+1);
t = linspace(0,t_f,t_n);

[X,Y] = meshgrid(x, y);

% loop to set initial conditions
for i = 1:n*L
    for j = 1:n*H
        u_initial(j,i) = x(i).*y(j).*(L-x(i)).*(H-y(j));
    end
end

for i = 1:t_n
    results_analytical(:,:,i) = u_a(X, Y, t(i));
end

figure()
% surf(u_initial) % display initial conditions for verification

% loop to find solutions, use if statements for boundary conditions
u_p = u_initial;
u_c = u_initial;
for k = 1:t_n
    for i = 1:n*L     
        for j = 1:n*H
            if i == 1 || j == 1 || i == n*L || j == n*H %% enforcing boundary conditions
                u_n(j,i) = 0;
            else
                u_n(j,i) = 2*(1-2*lambda2)*u_c(j,i) + lambda2*(u_c(j+1,i) + u_c(j-1,i) + u_c(j,i+1) + u_c(j,i-1)) - u_p(j,i);
            end         
        end
    end
    u_p = u_c;
    u_c = u_n;
    results_numerical(:,:,k) = u_n;
end

figure()

if animated_plot
    for time = 1:t_n
        subplot(1,2,1)
        surf(results_numerical(:,:,time));
        axis([0 n*L 0 n*H -2.5 2.5]);
        title('numerical solution')
        subplot(1,2,2)
        surf(results_analytical(:,:,time));
        axis([0 n*L 0 n*H -2.5 2.5]);
        title('analytical solution')
        pause(1/60)
    end
end

figure()
colour = 0;
num_plots = 10;
for time = 1:t_n/num_plots:t_n
    formatspec = 't = %.2f s';
    legendstr = sprintf(formatspec, time);
    plot(results_numerical(round(n*H/2), :, time), 'Color',[1-colour,0,colour], 'DisplayName', legendstr);
    colour = colour + 1/num_plots;
    hold on
    legend
end
title("u(x, H/2, t) for various times, numerical solution")

figure()

colour = 0;
num_plots = 10;
for time = 1:t_n/num_plots:t_n
    formatspec = 't = %.2f s';
    legendstr = sprintf(formatspec, time);
    plot(results_analytical(round(n*H/2), :, time), 'Color',[1-colour,0,colour], 'DisplayName', legendstr);
    colour = colour + 1/num_plots;
    hold on
    legend
end

title("u(x, H/2, t) for various times, analytical solution")