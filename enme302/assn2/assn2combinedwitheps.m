clc; close all; clear

% VARIABLES
L = 2;
H = 3;

g = 0;
s = 0;
c = 6;

animated_plot = 0;
plot_at_all = 0;
lambda = 0.6 % constant courant number
lambda2 = lambda^2;
deltas = linspace(0.5, 0.005, 50)
epss = []

for delta = deltas
grid_size = 1; % physical size of grid
t_f = 0.5; % end time of sim

n = round(grid_size/delta); % calculating spatial grid points to maintain constant courant number

dt = lambda*delta/c;
t_n = t_f/dt;

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
u_initial = [];
for i = 1:length(x)
    for j = 1:length(y)
        u_initial(j,i) = x(i).*y(j).*(L-x(i)).*(H-y(j));
    end
end

results_analytical = [];
for i = 1:t_n
    results_analytical(:,:,i) = u_a(X, Y, t(i));
end

% figure()
% surf(u_initial) % display initial conditions for verification

% loop to find solutions, use if statements for boundary conditions
u_p = u_initial;
u_c = u_initial;
results_numerical = [];
for k = 1:t_n
    for i = 1:length(x)     
        for j = 1:length(y)
            if i == 1 || j == 1 || i == length(x)  || j == length(y)  %% enforcing boundary conditions
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


if plot_at_all
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
colour = 0
num_plots = 10
for time = 1:t_n/num_plots:t_n
    formatspec = 't = %.2f s'
    legendstr = sprintf(formatspec, time)
    plot(results_numerical(round(n*H/2), :, time), 'Color',[1-colour,0,colour], 'DisplayName', legendstr);
    colour = colour + 1/num_plots
    hold on
    legend
end
title("u(x, H/2, t) for various times, numerical solution")

figure()

colour = 0
num_plots = 10
for time = 1:t_n/num_plots:t_n
    formatspec = 't = %.2f s'
    legendstr = sprintf(formatspec, time)
    plot(results_numerical(round(n*H/2), :, time), 'Color',[1-colour,0,colour], 'DisplayName', legendstr);
    colour = colour + 1/num_plots
    hold on
    legend
end

title("u(x, H/2, t) for various times, analytical solution")
end

% looping to calculate epsilon
sum = 0;
for i = 1:length(x)   
    for j = 1:length(y)
        sum = sum + ((results_numerical(j, i, end)) - results_analytical(j, i, end))^2;       
    end
end
eps = sqrt(delta^2*sum);
epss = [epss eps];
end
plot(deltas, epss)