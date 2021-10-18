clc; close all; clear

% VARIABLES
L = 2;
H = 3;

u_e = @(x, y, t) x.*(L-x).*y.*(H-y).*(1+1/2.*t);

c = 6;
g = 0;
s = 0;

lambda = 0.6; % constant courant number
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
t = linspace(0,t_f,t_n);

[X,Y] = meshgrid(x, y);

u_initial = u_e(X, Y, 0);

for i = 1:t_n
    results_analytical(:,:,i) = u_e(X, Y, t(i));
    surf(results_analytical(:,:,i));
end


% loop to find solutions, use if statements for boundary conditions
u_c = u_initial;
for k = 1:t_n
    for i = 1:length(x)-1    
        for j = 1:length(y)-1
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

figure()

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

figure()
results = []
for time = 1:t_n
    results = [results results_numerical(round(n*H/2), round(n*L/2), time)];
    clc
    fprintf("iteration %d of %d\n", time, t_n)
    fprintf("progress: %.2f%% \n", time/t_n*100)
end

plot(t, results)
