clc; close all; clear

% VARIABLES
L = 2;
H = 3;

u_e = @(x, y, t) x.*(L-x).*y.*(H-y).*(1+1/2.*t);

c = 6;
g = @(x, y, t) 0.5 .* u_e(x, y, t);
s = @(x, y, t) 2.*c^2 .* (1 + 0.5.*t) .* (y .* (H-y) + x .* (L-x)); 

lambda = 0.4; % constant courant number
lambda2 = lambda^2;

delta = 0.2;

animated_plot = 0;

grid_size = 1; % physical size of grid
t_f = 5; % end time of sim

n = round(grid_size/delta); % calculating spatial grid points to maintain constant courant number

dt = lambda*delta/c;
t_n = t_f/dt;

% ANALYTICAL SOLUTION

x = linspace(0,L,L/delta+1);
y = linspace(0,H,H/delta+1);
t = linspace(0,t_f,t_n);

[X,Y] = meshgrid(x, y);

u_initial = u_e(X, Y, 0);
results_analytical = zeros(length(y), length(x), length(t));
for i = 1:t_n
    results_analytical(:,:,i) = u_e(X, Y, t(i));
    surf(results_analytical(:,:,i));
end


% loop to find solutions, use if statements for boundary conditions
u_c = u_initial;
results_numerical = zeros(length(y), length(x), length(t));
for k = 1:t_n
    u_n = zeros(length(y), length(x));
    for i = 1:length(x)     
        for j = 1:length(y)
            if i == 1 || j == 1 || i == length(x) || j == length(y) %% enforcing boundary conditions
                u_n(j,i) = 0;
            else
                if k == 1 % implementing special case for first time step, using ghost node to implement V.N. B.C.
                    u_n(j,i) = 2*(1-2*lambda2)*u_c(j,i) + lambda2*(u_c(j+1,i) + u_c(j-1,i) + u_c(j,i+1) + u_c(j,i-1)) - (u_c(j,i)-2.*g(x(i), y(j), t(k))*dt) + dt^2 * s(x(i), y(j), t(k));
                else
                    u_n(j,i) = 2*(1-2*lambda2)*u_c(j,i) + lambda2*(u_c(j+1,i) + u_c(j-1,i) + u_c(j,i+1) + u_c(j,i-1)) - u_p(j,i) + dt^2 * s(x(i), y(j), t(k));
                end
            end         
        end
    end
    u_p = u_c;
    u_c = u_n;
    results_numerical(:,:,k) = u_n;
end

if animated_plot
filename = 'manufactured solutions.gif';
gif(filename)
for time = 1:t_n
    subplot(1,2,1)
    surf(X, Y, results_numerical(:,:,time));
    axis([0 L 0 H]);
    title('numerical solution')
    subplot(1,2,2)
    surf(X, Y, results_analytical(:,:,time));
    axis([0 L 0 H]);
    title('analytical solution')
    if time == 1
        gif(filename);
    else
        gif('DelayTime', 1/60)
    end
end
end

figure()
epsilon_array = zeros(1, length(t));
for time = 1:length(t)
    tosum = (results_numerical(:,:,time) - results_analytical(:,:,time)).^2;
    sum2 = sum(tosum, 'all');
    eps = sqrt(delta^2.*sum2);
    epsilon_array(time) = eps;
end

plot(t, epsilon_array)
xlabel('Time elapsed [s]')
ylabel('Error (epsilon) [m]')
title('Error against elapsed time')

