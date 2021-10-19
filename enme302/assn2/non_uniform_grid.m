clear; close all; clc

f = @(x, y) 5*normpdf(x,0,1);        % function for determining initial conditions
c = @(x, y) 5;       % function for determining variable wave speed
t_final = 2;            % [s], length over which simulation is run
time_points  = 1000;     % number of time points to iterate over
n = 10;                 % number of grid points


delta_t = t_final/time_points;  % [s], resulting time step for simulation

L = 2;                          % [m], x dimension of simulation space
H = 3;                          % [m], y dimension of simulation space

% PLACEHOLDER PENDING PROPER DISCRETISATION CODE
x = linspace(0,L,n*L);                      % vector of grid points in x
y = linspace(0,H,n*H);                      % vector of grid points in y
t = linspace(0,t_final,time_points);        % vector of time points


% NUMERICAL DISCRETISED SOLUTION
u = zeros(length(y), length(x), length(t));     % preallocate u for speed
for k = 1:time_points
    for j = 1:length(y)
        for i = 1:length(x)
            % enforcing boundary conditions
            if i == 1 || j == 1 || i == length(x) || j == length(y) 
                u(j,i,k) = 0;
            else
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
end
figure(1)
surf(u(:,:,1))

figure(2)
for i = 1:length(t)
    % s = surf(u(:,:,i), 'FaceAlpha', 1, 'EdgeColor', 'interp', 'FaceLighting', 'flat');
    s = surf(u(:,:,i), 'FaceAlpha', 1, 'FaceLighting', 'flat');
    axis([0 n*L 0 n*H -60 60]);
    pause(1/144)
end
