clc; close all; clear

% VARIABLES
L = 2;
H = 3;

u_e = @(x, y, t) x.*(L-x).*y.*(H-y).*(1+1/2.*t);

c = 6;
g = @(x, y, t) 0.5 .* u_e(x, y, t);
s = @(x, y, t) 2.*c^2 .* (1 + 0.5.*t) .* (y .* (H-y) + x .* (L-x)); 

lambda = 0.25; % constant courant number
lambda2 = lambda^2;

delta = 0.1;

animated_plot = 0;

grid_size = 1; % physical size of grid
t_f = 1; % end time of sim

n = round(grid_size/delta); % calculating spatial grid points to maintain constant courant number

dt = lambda*delta/c;
t_n = t_f/dt;

% ANALYTICAL SOLUTION

x = linspace(0,L,L/delta+1);
y = linspace(0,H,H/delta+1);
t = 0:dt:t_f;

[X,Y] = meshgrid(x, y);

u_initial = u_e(X, Y, 0);
results_analytical = zeros(length(y), length(x), length(t));
for i = 1:length(t)
    results_analytical(:,:,i) = u_e(X, Y, t(i));
end

surf(X,Y,results_analytical(:,:,1))

% loop to find solutions, use if statements for boundary conditions
u = zeros(length(y), length(x), length(t));     % preallocate u for speed
u(:,:,1) = u_initial;
u(:,:,2) = u_initial;
for k = 2:length(t)
    for j = 1:length(y)
        for i = 1:length(x)
            % enforcing boundary conditions
            if i == 1 || j == 1 || i == length(x) || j == length(y) 
                u(j,i,k+1) = 0;
            else
                if k == 1
                    % THIS LINE SUCKS AND IT'S WRONG and i hate it
                    u(j,i,k+1) = (1-2*lambda2)*u(j,i,k)+0.5*lambda2*(u(j,i+1,k)+u(j,i-1,k)+u(j+1,i,k)+u(j-1,i,k));
                else
                    u(j,i,k+1) = 2*(1-2*lambda2)*u(j,i,k)+lambda2*(u(j,i+1,k)+u(j,i-1,k)+u(j+1,i,k)+u(j-1,i,k))-u(j,i,k-1);
                end
            end
        end
    end
end

results_numerical = u;

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

