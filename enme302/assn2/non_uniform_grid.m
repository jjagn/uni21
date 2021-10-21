clear; close all; clc

f = @(x, y) sin(x/2)+0.1;        % function for determining initial conditions
c = @(x, y) x^2 + 1;       % function for determining variable wave speed
t_final = 5;            % [s], length over which simulation is run
time_points = 1000;     % number of time points to iterate over
n = 10;                 % number of grid points

delta_t = t_final/time_points;  % [s], resulting time step for simulation

L = 2;                          % [m], x dimension of simulation space
H = 3;                          % [m], y dimension of simulation space
alpha = 0.3;

% PLACEHOLDER PENDING PROPER DISCRETISATION CODE
x = generate_grid(L, n*L, alpha);           % vector of grid points in x
y = generate_grid(H, n*H, alpha);           % vector of grid points in y
t = linspace(0,t_final,time_points);        % vector of time points

[X,Y] = meshgrid(x, y);


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
surf(X, Y, u(:,:,1))

figure(2)
for i = 1:length(t)
    % s = surf(u(:,:,i), 'FaceAlpha', 1, 'EdgeColor', 'interp', 'FaceLighting', 'flat');
    s = surf(X,Y,u(:,:,i), 'FaceAlpha', 1, 'FaceLighting', 'flat');
    axis([0 L 0 H]);
    pause(1/144)
end

function grid_points = generate_grid(L,N,alpha)
x_k = @(alpha,K,L,n) L*alpha*sum((1+alpha).^(K-1))/(2*((1+alpha)^n-1));

n = N/2;
X = zeros(1,n);
for i  = 1:n-1
    K = 1:i;
    X(i+1) = x_k(alpha, K, L, n);
end

grid_points = zeros(1,N);
grid_points(1, 1:n) = X;
grid_points(1, n+1:N) = L-fliplr(X);
end


