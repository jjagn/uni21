clc; close all; clear

% VARIABLES
L = 2;
H = 3;

g = 0;
s = 0;
c = 6;

lambda = 0.25; % constant courant number
lambda2 = lambda^2;
nx = [2 4 8 16 32 64 128];
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


x = linspace(0,L,n*L);
y = linspace(0,H,n*H);
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
figure()
surf(u_initial)

results_analytical = zeros(length(y), length(x), length(t));
for i = 1:t_n
    results_analytical(:,:,i) = u_a(X, Y, t(i));
end

% loop to find solutions, use if statements for boundary conditions
u_p = u_initial;
u_c = u_initial;
results_numerical = zeros(length(y), length(x), length(t));
for k = 1:t_n
    u_n = zeros(length(y), length(x));
    for i = 1:length(x)     
        for j = 1:length(y)
            if i == 1 || j == 1 || i == length(x)  || j == length(y)  %% enforcing boundary conditions
                u_n(j,i) = 0;
            else
                u_n(j,i) = 2*(1-2*lambda2)*u_c(j,i) + lambda2*(u_c(j+1,i) + u_c(j-1,i) + u_c(j,i+1) + u_c(j,i-1)) - u_p(j,i);
            end         
        end
    end
    u_p = u_c;
    u_c = u_n;
    results_numerical(:,:,k) = u_n;
end

% calculating epsilon
tosum = (results_numerical(:,:,end) - results_analytical(:,:,end)).^2;
sum2 = sum(tosum, 'all');

eps = sqrt(delta^2.*sum2);
epss = [epss eps];
end
plot(fliplr(deltas), epss)