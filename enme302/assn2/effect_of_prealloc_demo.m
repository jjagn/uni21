clear; close all; clc

f = @(x, y) x-y;        % function for determining initial conditions
c = @(x, y) 2.*x;       % function for determining variable wave speed
t_final = 1;            % [s], length over which simulation is run
time_points  = 100;     % number of time points to iterate over
n = 10;                 % number of grid points


delta_t = t_final/time_points;  % [s], resulting time step for simulation

L = 2;                          % [m], x dimension of simulation space
H = 3;                          % [m], y dimension of simulation space

% PLACEHOLDER PENDING PROPER DISCRETISATION CODE
x = linspace(0,L,n*L);                      % vector of grid points in x
y = linspace(0,H,n*H);                      % vector of grid points in y
t = linspace(0,t_final,time_points);        % vector of time points

elapsed = []
for iter = 1:100
% NUMERICAL DISCRETISED SOLUTION
tic
%u = zeros(length(y), length(x), length(t));     % preallocate u for speed
for k = 1:t_final
    for j = 1:length(y)
        for i = 1:length(x)
            % enforcing boundary conditions
            if i == 1 || j == 1 || i == length(x) || j == length(y) 
                u(j,i,k) = 0;
            end
            % setting initial conditions, CURRENTLY IGNORING THE EFFECTS OF
            % VON NEUMANN BOUNDARY CONDITION G
            if k == 1 
                u(j,i,k) = f(x(i), y(j));
            else
            u(j,i,k+1) = c(x(i),y(i))^2*delta_t^2*...
                (((u(j,i+1,k)-2*u(j,i,k)+u(j,i-1,k))/(x(i+1)-x(i))^2) + ...
                ((u(j+1,i,k)-2*u(j,i,k)+u(j-1,i,k))/(y(i+1)-y(i))^2)) + ...
                2*u(j,i,k)-u(j,i,k-1);
            end
        end
    end
end
elapsed = [elapsed toc];
end
plot(elapsed)
mean = mean(elapsed);
min = min(elapsed);
max = max(elapsed);
median = median(elapsed);
fprintf("mean: %.5f\n", mean)
fprintf("min: %.5f\n", min)
fprintf("max: %.5f\n", max)
fprintf("median: %.5f\n", median)

