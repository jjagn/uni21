clear; close all; clc

L = 2;                          % [m], x dimension of simulation space
H = 3;                          % [m], y dimension of simulation space
alpha = 0.05;

f = @(x, y) sin(5*x)+ cos(y-(H/2));        % function for determining initial conditions
c = @(x, y) sin(x)+2;       % function for determining variable wave speed
t_final = 1;            % [s], length over which simulation is run
time_points = 1000;     % number of time points to iterate over
n = 30;                 % number of grid points

delta_t = t_final/time_points;  % [s], resulting time step for simulation



% DISCRETISATION CODE FOR NON-UNIFORM GRID
x = generate_grid(L, n*L, alpha);           % vector of grid points in x
y = generate_grid(H, n*H, alpha);           % vector of grid points in y
t = linspace(0,t_final,time_points);        % vector of time points

[X,Y] = meshgrid(x, y);

lambda = zeros(1,length(x));
for i = 1:length(x)-1
    for j = 1:length(y)
        lambda(i) = delta_t*c(x(i),y(j))/(x(i+1)-x(i));
    end
end

% NUMERICAL DISCRETISED SOLUTION
u = zeros(length(y), length(x), length(t));     % preallocate u for speed
u(:,:,1) = f(X,Y); 
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
                    dx = x(i+1) - x(i);
                    dy = y(j+1) - y(j);
        
                    c_iph = (c(x(i+1),y(i))^2+c(x(i),y(i))^2)/2;
                    c_imh = (c(x(i-1),y(j))^2+c(x(i),y(j))^2)/2;
                    c_jph = (c(x(i),y(j+1))^2+c(x(i),y(j))^2)/2;
                    c_jmh = (c(x(i),y(j-1))^2+c(x(i),y(j))^2)/2;
                    
                    denominator_i = c_iph*u(j,i+1,k)-(c_iph+c_imh)*u(j,i,k)+c_imh*u(j,i-1,k);
                    
                    denominator_j = c_jph*u(j+1,i,k)-(c_jph+c_jmh)*u(j,i,k)+c_jmh*u(j-1,i,k);
                    
                    inside = denominator_i/dx^2 + denominator_j/dy^2;
                    
                    u(j,i,k+1) = 0.5*(delta_t^2*inside)+u(j,i,k);
                else
                
                    dx = x(i+1) - x(i-1);
                    dy = y(j+1) - y(j-1);
        
                    c_iph = (c(x(i+1),y(i))^2+c(x(i),y(i))^2)/2;
                    c_imh = (c(x(i-1),y(j))^2+c(x(i),y(j))^2)/2;
                    c_jph = (c(x(i),y(j+1))^2+c(x(i),y(j))^2)/2;
                    c_jmh = (c(x(i),y(j-1))^2+c(x(i),y(j))^2)/2;
                    
                    denominator_i = c_iph*u(j,i+1,k)-(c_iph+c_imh)*u(j,i,k)+c_imh*u(j,i-1,k);
                    
                    denominator_j = c_jph*u(j+1,i,k)-(c_jph+c_jmh)*u(j,i,k)+c_jmh*u(j-1,i,k);
                    
                    inside = denominator_i/dx^2 + denominator_j/dy^2;
                    
                    u(j,i,k+1) = delta_t^2*inside+2*u(j,i,k)-u(j,i,k-1);
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


