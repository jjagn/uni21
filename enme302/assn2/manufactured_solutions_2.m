%% RESET
clear; close all; clc

%% GENERAL SETUP
% METHOD OF MANUFACTURED SOLUTIONS V2
% VARIABLES

% VARIABLES
L = 2;
H = 3;

u_e = @(x, y, t) x.*(L-x).*y.*(H-y).*(1+0.5.*t);

c = 6;

f = @(x,y) u_e(x,y,0);
g = @(x,y) 0.5.* x.*(L-x).*y.*(H-y);

s = @(x,y,t) 2.*c^2 .* (1+0.5.*t) .* (y.*(H-y)+x.*(L-x)); 

lambda = 0.5; % constant courant number
lambda2 = lambda^2;

delta = 0.1;

animated_plot = 0;

grid_size = 1; % physical size of grid
t_f = 100; % end time of sim

n = round(grid_size/delta); % calculating spatial grid points to maintain constant courant number

dt = lambda*delta/c;
t_n = t_f/dt;

%% SOLUTION SETUP
x = 0:delta:L;
y = 0:delta:H;
t = 0:dt:t_f;

[X,Y] = meshgrid(x,y);
x_len = length(x);
y_len = length(y);
t_len = length(t);

%% ANALYTICAL SOLVER
results_analytical = zeros(y_len, x_len, t_len);
for i = 1:t_len
    results_analytical(:,:,i) = u_e(X, Y, t(i));
end
%% NUMERICAL SOLUTION
u = zeros(y_len,x_len,t_len); % prealloc for speed
u(:,:,1) = f(X,Y); % set initial conditions

for k = 1:t_len-1
    for j = 1:y_len
        for i = 1:x_len
            % capture boundary conditions
            if i == 1 || i == x_len || j == 1 || j == y_len
                u(j,i,k) = 0;
            elseif k == 1 % special case for first time step
                u(j,i,k+1) = (1-2*lambda2)*u(j,i,k)+0.5*lambda2*...
                    (u(j,i+1,k)+u(j+1,i,k)+u(j,i-1,k)+u(j-1,i,k))+...
                    0.5*s(x(i),y(j),t(k))*dt^2+dt*g(x(i),y(j));
            else % normal operation for internal nodes after first step
                u(j,i,k+1) = 2*(1-2*lambda2)*u(j,i,k)+lambda2*...
                    (u(j,i+1,k)+u(j+1,i,k)+u(j,i-1,k)+u(j-1,i,k))+...
                    s(x(i),y(j),t(k))*dt^2-u(j,i,k-1);
            end
        end
    end
end

%% ANIMATED PLOT FOR NUMERICAL SOLUTION
figure()
for i = 1:t_len
    surf(X,Y,u(:,:,i))
    pause(1/60)
end  

%% ANIMATED SURFACE OF ANALYTICAL SOLUTION
figure()
for i = 1:t_len
    surf(X,Y,results_analytical(:,:,i))
    pause(1/60)
end

%% ANIMATED PLOT FOR BOTH
figure()
for i = 1:t_len
    subplot(1,2,1)
    surf(X,Y,u(:,:,i))
    title('numerical solution')
    axis([0 L 0 H 0 4])
    subplot(1,2,2)
    surf(X,Y,results_analytical(:,:,i))
    axis([0 L 0 H 0 4])
    title('analytical solution')
    pause(1/60)
end  

%% FUCK IT, HAVE A GIF
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
filename = 'q3animated.gif';

maxmax = max(max(max(u)));
minmin = min(min(min(u)));

for i = 1:length(t)
    if i == 1
        gif(filename)
    else
        gif('DelayTime', 1/60)
    end
    subplot(4,4,[1:2,5:6])
    surf(X,Y,u(:,:,i))
    title('Numerical solution')
    axis([0 L 0 H minmin maxmax]);
    caxis([minmin maxmax])
    colorbar
    subplot(4,4,[3:4,7:8])
    surf(X,Y,results_analytical(:,:,i))
    axis([0 L 0 H minmin maxmax]);
    caxis([minmin maxmax])
    colorbar
    title('Analytical solution')
    subplot(4,4,[10:11,14:15])
    error = abs(u(:,:,i)-results_analytical(:,:,i));
    surf(X,Y,error);
    axis([0 L 0 H]);
    colorbar
    title('Error')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end 

%% ERROR CALCULATION
epsilon_array = zeros(1, length(t));
for time = 1:length(t)
    tosum = (u(:,:,time) - results_analytical(:,:,time)).^2;
    sum2 = sum(tosum, 'all');
    eps = sqrt(delta^2.*sum2);
    epsilon_array(time) = eps;
end

plot(t, epsilon_array)
xlabel('Time elapsed [s]')
ylabel('Error (epsilon) [m]')
title('Error against elapsed time')
