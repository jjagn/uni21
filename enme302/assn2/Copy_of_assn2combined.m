clc; close all; clear

% VARIABLES
L = 2;
H = 3;

g = 0;
s = 0;
c = 6;
n = 10; % spatial grid points
grid_size = 1; % physical size of grid
t_n = 200; % number of time points
t_f = 1; % end time of sim

animated_plot = 1;

delta = grid_size/n;
dt = t_f/t_n;

lambda = c*dt/delta;
lambda2 = lambda*lambda;

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
u_initial = zeros(length(y), length(x));
for i = 1:n*L
    for j = 1:n*H
        u_initial(j,i) = x(i).*y(j).*(L-x(i)).*(H-y(j));
    end
end

results_analytical = zeros(length(y), length(x), length(t));
for i = 1:t_n
    results_analytical(:,:,i) = u_a(X, Y, t(i));
end

figure()
%surf(X,Y, u_initial) % display initial conditions for verification

% loop to find solutions, use if statements for boundary conditions
results_numerical = zeros(length(y), length(x), length(t));
results_numerical(:,:,1) = u_initial;
surf(x,Y,results_numerical(:,:,1))
for k = 1:t_n
    for i = 1:length(x)-1     
        for j = length(y)-1
            if i == 1 || j == 1 || i == length(x) || j == length(y) %% enforcing boundary conditions
                results_numerical(j,i,k+1) = 0
            elseif k == 1
                results_numerical(j,i,k+1) = 2*(1-2*lambda2).*results_numerical(j,i,k) + lambda2.*(results_numerical(j+1,i,k) + results_numerical(j-1,i,k) + results_numerical(j,i+1,k) + results_numerical(j,i-1,k));
            else
                results_numerical(j,i,k+1) = 2*(1-2*lambda2).*results_numerical(j,i,k) + lambda2.*(results_numerical(j+1,i,k) + results_numerical(j-1,i,k) + results_numerical(j,i+1,k) + results_numerical(j,i-1,k)) - results_numerical(j,i,k-1);
            end         
        end
    end
end

figure()

if animated_plot
    for time = 1:t_n
        subplot(1,2,1)
        surf(X, Y, results_numerical(:,:,time));
        axis([0 L 0 H -2.5 2.5]);
        title('numerical solution')
        subplot(1,2,2)
        surf(X, Y, results_analytical(:,:,time));
        axis([0 L 0 H -2.5 2.5]);
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