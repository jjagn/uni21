clc; close all; clear

% VARIABLES
L = 2;
H = 3;

g = 0;
s = 0;
c = 6;

lambda = 0.25; % constant courant number
lambda2 = lambda^2;
nx = [8 12 16 24 32 48 56 64 72 96 108 128 192 256];
nx = fliplr(nx)
% n_iterations = 20
% nx = linspace(4, 4*n_iterations, n_iterations)

epss = [];
deltas = [];
for n = nx
grid_size = 1; % physical size of grid
t_f = 0.5; % end time of sim

ny = H/L.* n;
delta = L/n;
deltas = [deltas delta];
dt = lambda*delta/c;
t_n = t_f/dt;

% ANALYTICAL SOLUTION
taylor_expansions = 5;

u_a = @(x, y, t) 0;

for b = 1:taylor_expansions
    for a = 1:taylor_expansions
        u_a = @(x, y, t) u_a(x, y, t) + (576/pi^6) .* (((1 + (-1)^(a+1)) * (1+(-1)^(b+1))) / (a^3*b^3)).*sin(a.*pi./2.*x).*sin(b.*pi./3.*y).*cos(pi.*sqrt(9.*a.^2+4.*b.^2).*t);
    end
end


x = linspace(0,L,n);
y = linspace(0,H,ny);
t = linspace(0,t_f,t_n);

[X,Y] = meshgrid(x, y);

% loop to set initial conditions
u_initial = zeros(length(y), length(x));
length(x)
length(y)
for i = 1:length(x)
    for j = 1:length(y)
        u_initial(j,i) = x(i).*y(j).*(L-x(i)).*(H-y(j));
    end
end
% figure()
% surf(X,Y,u_initial)

results_analytical = zeros(length(y), length(x), length(t));
for i = 1:t_n
    results_analytical(:,:,i) = u_a(X, Y, t(i));
end

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
                u(j,i,k+1) = 2*(1-2*lambda2)*u(j,i,k)+lambda2*(u(j,i+1,k)+u(j,i-1,k)+u(j+1,i,k)+u(j-1,i,k))-u(j,i,k-1);
            end
        end
    end
end

% calculating epsilon
% figure()
% surf(X,Y,u(:,:,end))
% title('numerical')
% figure()
% surf(X,Y,results_analytical(:,:,end))
% title('analytical')
tosum = (u(:,:,end-1) - results_analytical(:,:,end)).^2;
sum2 = sum(tosum, 'all');

eps = sqrt(delta^2.*sum2);
epss = [epss eps];
end
plot(deltas, epss)
title('error against grid spacing')
xlabel('grid spacing (delta)')
ylabel('error (epsilon)')