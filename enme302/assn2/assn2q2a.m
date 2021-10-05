clc; close all; clear

% VARIABLES
L = 2;
H = 3;

g = 0;
s = 0;
c = 6;
n = 50; % spatial grid points
grid_size = 1; % physical size of grid
t_n = 500; % number of time points
t_f = 1; % end time of sim

delta = grid_size/n;
dt = t_f/t_n;

lambda = c*dt/delta
lambda2 = lambda*lambda


x = linspace(0,L,L/delta+1);
y = linspace(0,H,H/delta+1);
t = linspace(0,t_f,t_n);

% loop to set initial conditions
for i = 1:n*L
    for j = 1:n*H
        u_initial(i,j) = x(i).*y(j).*(L-x(i)).*(H-y(j));
    end
end

figure()
surf(u_initial) % display initial conditions for verification

% loop to find solutions, use if statements for boundary conditions
u_p = u_initial;
u_c = u_initial;
for k = 1:t_n
    for i = 1:n*L     
        for j = 1:n*H
            if i == 1 || j == 1 || i == n*L || j == n*H %% enforcing boundary conditions
                u_n(i,j) = 0;
            else
                u_n(i,j) = 2*(1-2*lambda2)*u_c(i,j) + lambda2*(u_c(i+1,j) + u_c(i-1,j) + u_c(i,j+1) + u_c(i,j-1)) - u_p(i,j);
            end         
        end
    end
    u_p = u_c;
    u_c = u_n;
    solution(:,:,k) = u_n;
end

filename = 'numerical_soln.gif';
gif(filename)
for time = 1:t_n
    surf(solution(:,:,time));
    axis([0 n*H 0 n*L -2.5 2.5]);
    % pause(1/60)
    if time == 1
        gif(filename)
    else 
        gif('DelayTime', 1/60)
    end
end