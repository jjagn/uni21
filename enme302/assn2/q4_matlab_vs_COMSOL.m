clc; close all; clear

% VARIABLES
L = 2;
H = 3;

u_e = @(x, y, t) x.*(L-x).*y.*(H-y).*(1+1/2.*t);

c = 6;
g = @(x, y, t) 0;
s = @(x, y, t) 0;

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
a_lim = 5;      % number of fourier series terms
b_lim = 5;      % "" ""

u_a = @(x, y, t) 0;

for b = 1:1:b_lim
    for a = 1:1:a_lim
        u_a = @(x, y, t) u_a(x, y, t) + (576/pi^6) .* (((1 + (-1)^(a+1)) * (1+(-1)^(b+1))) / (a^3*b^3)).*sin(a.*pi./2.*x).*sin(b.*pi./3.*y).*cos(pi.*sqrt(9.*a.^2+4.*b.^2).*t);
    end
end

x = linspace(0,L,L/delta+1);
y = linspace(0,H,H/delta+1);
t = linspace(0,t_f,t_n);

[X,Y] = meshgrid(x, y);

% loop to set initial conditions
for i = 1:length(x)
    for j = 1:length(y)
        u_initial(j,i) = x(i).*y(j).*(L-x(i)).*(H-y(j));
    end
end

for i = 1:t_n
    results_analytical(:,:,i) = u_a(X, Y, t(i));
end


% loop to find solutions, use if statements for boundary conditions
u_c = u_initial;
for k = 1:t_n
    for i = 1:length(x)   
        for j = 1:length(y)
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

% for time = 1:t_n
%     subplot(1,2,1)
%     surf(results_numerical(:,:,time));
%     axis([0 n*L 0 n*H -2.5 2.5]);
%     title('numerical solution')
%     subplot(1,2,2)
%     surf(results_analytical(:,:,time));
%     axis([0 n*L 0 n*H -2.5 2.5]);
%     title('analytical solution')
%     pause(1/60)
% end

d=readmatrix('midpoint.txt','headerlines',9, 'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join');
d = d(4:end)

figure()
results = []
for time = 1:t_n
    results = [results results_numerical(round(n*H/2), round(n*L/2), time)];
    clc
    fprintf("iteration %d of %d\n", time, t_n)
    fprintf("progress: %.2f%% \n", time/t_n*100)
end

plot(t, results)
hold on
plot(t, d)
legend('MATLAB', 'COMSOL')
title('u(L/2, H/2, t) for MATLAB and COMSOL numerical solutions')
xlabel('Time [s]'); ylabel('Displacement (u) [m]')