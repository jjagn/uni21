%% RESET
clear; close all; clc

%% ERROR AGAINST GRID SPACING
% set true if trying to collect error against grid spacing
iterate_over_multiple_deltas = false

if iterate_over_multiple_deltas
    % create array of deltas
    deltas = [0.01 0.02 0.025 0.05 0.1 0.2 0.25 0.5 1]; 
    epsilon_array = [];
else
    deltas = 0.05; % sets grid spacing
end

for delta = deltas
%% GENERAL SETUP
% METHOD OF MANUFACTURED SOLUTIONS V2
% VARIABLES

% VARIABLES
L = 2;      % domain length [m]
H = 3;      

c = 6;

f = @(x,y) x.*y.*(L-x).*(H-y);
g = @(x,y) 0;
s = @(x,y,t) 0;

lambda = 0.5; % constant courant number
lambda2 = lambda^2;

animated_plot = 0;

% calculate other params from courant number

grid_size = 1; % physical size of grid
t_f = 1; % end time of sim

n = round(grid_size/delta); % calculating spatial grid points to maintain constant courant number

dt = lambda*delta/c;
t_n = t_f/dt;

%% SOLUTION SETUP
 
% create arrays representing domain
x = 0:delta:L;
y = 0:delta:H;
t = 0:dt:t_f;

% meshgrid for vectorising
[X,Y] = meshgrid(x,y);
x_len = length(x);
y_len = length(y);
t_len = length(t);

%% ANALYTICAL SOLVER

% 5 taylor expansions seems to be enough
taylor_expansions = 5;

u_a = @(x, y, t) 0;

for b = 1:taylor_expansions
    for a = 1:taylor_expansions
        % update function handle each loop
        u_a = @(x, y, t) u_a(x, y, t) + (576/pi^6) .* (((1 + (-1)^(a+1)) * (1+(-1)^(b+1))) / (a^3*b^3)).*sin(a.*pi./2.*x).*sin(b.*pi./3.*y).*cos(pi.*sqrt(9.*a.^2+4.*b.^2).*t);
    end
end
% prealloc for speed
results_analytical = zeros(y_len, x_len, t_len);
for i = 1:t_len
    results_analytical(:,:,i) = u_a(X, Y, t(i));
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
                    s(x(i),y(j),t(k))*dt^2+dt*g(x(i),y(j));
            else % normal operation for internal nodes after first step
                u(j,i,k+1) = 2*(1-2*lambda2)*u(j,i,k)+lambda2*...
                    (u(j,i+1,k)+u(j+1,i,k)+u(j,i-1,k)+u(j-1,i,k))+...
                    s(x(i),y(j),t(k))*dt^2-u(j,i,k-1);
            end
        end
    end
end

if iterate_over_multiple_deltas
    % summing error
    tosum = (u(:,:,end) - results_analytical(:,:,end)).^2;
    sum2 = sum(tosum, 'all');
    eps = sqrt(delta^2.*sum2);
    epsilon_array = [epsilon_array eps];
end

end

%% ANIMATED PLOT FOR NUMERICAL SOLUTION
figure()
for i = 1:t_len
    % plot numerical solution over time
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
    subplot(1,3,1)
    surf(X,Y,u(:,:,i))
    title('Numerical solution')
    axis([0 L 0 H -2.5 2.5])
    subplot(1,3,2)
    surf(X,Y,results_analytical(:,:,i))
    axis([0 L 0 H -2.5 2.5])
    title('Analytical solution')
    subplot(1,3,3)
    error = abs(u(:,:,i)-results_analytical(:,:,i));
    surf(X,Y,error);
    axis([0 L 0 H 0 0.0025])
    title('Error')
    pause(1/60)
end  
%% ERROR CALCULATION
epsilon_array = zeros(1, length(t));
for time = 1:length(t)
    tosum = (u(:,:,time) - results_analytical(:,:,time)).^2;
    sum2 = sum(tosum, 'all');
    eps = sqrt(delta^2.*sum2);
    epsilon_array(time) = eps;
end

%% MAKE GIF
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
filename = 'numerical_analytical_error.gif';
for i = 1:t_len
    if i == 1
        gif(filename)
    else
        gif('DelayTime', 1/60)
    end
    subplot(4,4,[1:2,5:6])
    surf(X,Y,u(:,:,i))
    title('Numerical solution')
    axis([0 L 0 H -2.5 2.5])
    caxis([-2.5 2.5])
    colorbar
    subplot(4,4,[3:4,7:8])
    surf(X,Y,results_analytical(:,:,i))
    axis([0 L 0 H -2.5 2.5])
    caxis([-2.5 2.5])
    colorbar
    title('Analytical solution')
    subplot(4,4,[10:11,14:15])
    error = abs(u(:,:,i)-results_analytical(:,:,i));
    surf(X,Y,error);
    axis([0 L 0 H 0 0.025])
    caxis([0 0.025])
    colorbar
    title('Error')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end 

%% PLOT ERROR
plot(deltas, epsilon_array)
xlabel('Grid spacing [m]')
ylabel('Error (epsilon) [m]')
title('Error against grid spacing')
