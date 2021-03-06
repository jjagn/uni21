clc; close all; clear

% VARIABLES
L = 2;
H = 3;

g = 0;
s = 0;
c = 6;
n = 10; % spatial grid points
grid_size = 1; % physical size of grid
t_n = 250; % number of time points
t_f = 1; % end time of sim

animated_plot = 1;
gif = 0;

delta = grid_size/n;
dt = t_f/t_n;

lambda = c*dt/delta
lambda2 = lambda*lambda;

% ANALYTICAL SOLUTION
fourier_limit = 5;

u_a = @(x, y, t) 0;

for b = 1:fourier_limit
    for a = 1:fourier_limit
        u_a = @(x, y, t) u_a(x, y, t) + (576/pi^6) .* (((1 + (-1)^(a+1)) * (1+(-1)^(b+1))) / (a^3*b^3)).*sin(a.*pi./2.*x).*sin(b.*pi./3.*y).*cos(pi.*sqrt(9.*a.^2+4.*b.^2).*t);
    end
end

x = linspace(0,L,L/delta+1);
y = linspace(0,H,H/delta+1);
t = linspace(0,t_f,t_n);

[X,Y] = meshgrid(x, y);

% loop to set initial conditions
u_initial = zeros(length(y), length(x));
for i = 1:n*L
    for j = 1:n*H
        u_initial(j,i) = x(i).*y(j).*(L-x(i)).*(H-y(j));
    end
end

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

figure()

errors = zeros(length(y), length(x), length(t));
for time = 1:length(t)
    errors(:,:,time) = results_analytical(:,:,time) - u(:,:,time);
end

if animated_plot
    if gif
    filename = 'combined_solutions.gif';
    gif(filename)
    end
    for time = 1:t_n
        
        subplot(1,3,1)
        surf(X, Y, u(:,:,time));
        axis([0 L 0 H -2.5 2.5]);
        caxis([-2.5 2.5])
%         colorbar
        title('numerical solution')
        subplot(1,3,2)
        surf(X, Y, results_analytical(:,:,time));
        axis([0 L 0 H -2.5 2.5]);
        caxis([-2.5 2.5])
%         colorbar
        title('analytical solution')
        subplot(1,3,3)
        surf(X,Y,errors(:,:,time))
        axis([0 L 0 H -0.5 0.5]);
        caxis([-2.5 2.5])
%         colorbar
        title('error')
        set(gcf, 'Position', [0, 0, 720, 480])
        if gif
        if time == 1
            gif(filename)
        else 
            gif('DelayTime', 1/60)
        end
        else
            pause(1/240)
        end
        
    end
end

figure()
colour = 0;
num_plots = 10;
for time = 1:t_n/num_plots:t_n
    formatspec = 't = %.2f s';
    legendstr = sprintf(formatspec, time);
    subplot(1,2,1)
    plot(u(round(n*H/2), :, time), 'Color',[1-colour,0,colour], 'DisplayName', legendstr);
    hold on
    legend
    title("u(x, H/2, t) for various times, numerical solution")
    subplot(1,2,2)
    plot(results_analytical(round(n*H/2), :, time), 'Color',[1-colour,0,colour], 'DisplayName', legendstr);
    hold on
    legend
    title("u(x, H/2, t) for various times, analytical solution")
    colour = colour + 1/num_plots;
end