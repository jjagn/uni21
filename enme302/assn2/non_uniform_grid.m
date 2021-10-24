%% RESET
clear; close all; clc

%% VARIABLES SETUP
L = 2;                          % [m], x dimension of simulation space
H = 3;                          % [m], y dimension of simulation space
alpha = 0.05;                   % scaling variable for grid gen

% this one is kinda cool
f = @(x, y) 1-sin(x);        % function for determining initial conditions
c = @(x, y) 2*x+0.5;           % function for determining variable wave speed
t_final = 1;                    % [s], length over which simulation is run
time_points = 1000;             % number of time points to iterate over
n = 10;                         % number of grid points

delta_t = t_final/time_points;  % [s], resulting time step for simulation

animated_plot = 1;
do_gif = 1;

% 10^-4 tolerance for comsol

%% GENERATE GRID

% DISCRETISATION CODE FOR NON-UNIFORM GRID
x = generate_grid(L, n*L, alpha);           % vector of grid points in x
y = generate_grid(H, n*H, alpha);           % vector of grid points in y
t = 0:delta_t:t_final;        % vector of time points

[X,Y] = meshgrid(x, y);

%% PLOT MESH
% surf(X,Y,X*0)

%% CHECK LAMBDA
% CHECK LAMBDA SO SOLUTION REMAINS NUMERICALLY STABLE (DEPRECATED)
lambda = zeros(1,length(x));
for i = 1:length(x)-1
    for j = 1:length(y)
        lambda(i) = delta_t*c(x(i),y(j))/(x(i+1)-x(i));
    end
end
lambda_max = max(lambda)

%% NUMERICAL DISCRETISED SOLUTION
u = zeros(length(y), length(x), length(t));     % preallocate u for speed
u(:,:,1) = f(X,Y); 
for k = 1:time_points
    for j = 1:length(y)
        for i = 1:length(x)
            % enforcing boundary conditions
            if i == 1 || j == 1 || i == length(x) || j == length(y) 
                u(j,i,k) = 0;
            else
                % setting initial conditions
                if k == 1 
                    dx = x(i+1) - x(i);
                    dy = y(j+1) - y(j);
                    
                    % precalculate discretised values for c
                    c_iph = (c(x(i+1),y(i))^2+c(x(i),y(i))^2)/2;
                    c_imh = (c(x(i-1),y(j))^2+c(x(i),y(j))^2)/2;
                    c_jph = (c(x(i),y(j+1))^2+c(x(i),y(j))^2)/2;
                    c_jmh = (c(x(i),y(j-1))^2+c(x(i),y(j))^2)/2;
                    
                    % more simplifying so that the last line isn't a hot
                    % mess
                    denominator_i = c_iph*u(j,i+1,k)-(c_iph+c_imh)*u(j,i,k)+c_imh*u(j,i-1,k);
                    
                    denominator_j = c_jph*u(j+1,i,k)-(c_jph+c_jmh)*u(j,i,k)+c_jmh*u(j-1,i,k);
                    
                    inside = denominator_i/dx^2 + denominator_j/dy^2;
                    
                    u(j,i,k+1) = 0.5*(delta_t^2*inside)+u(j,i,k);
                else
                    % NORMAL OPERATION
                    dx = x(i+1) - x(i);
                    dy = y(j+1) - y(j);
        
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
% DISPLAY INITIAL CONDITIONS
figure(1)
surf(X, Y, u(:,:,1))

%% PLOTTING ANIMATED GIF AND ANIMATING SOLUTION DEPENDING ON SETTINGS
figure()
do_gif = false;
if do_gif
    filename = 'matlabfinalnotryans.gif';
end
maxmax = max(max(max(u)));
minmin = min(min(min(u)));
for i = 1:length(t)
    % s = surf(u(:,:,i), 'FaceAlpha', 1, 'EdgeColor', 'interp', 'FaceLighting', 'flat');
    s = surf(X,Y,u(:,:,i), 'FaceAlpha', 1, 'FaceLighting', 'flat');
    formatspec = 'Time = %f seconds';
    titlestring = sprintf(formatspec, t(i));
    title(titlestring)
    axis([0 L 0 H minmin maxmax]);
    caxis([minmin maxmax])
    colorbar
    if do_gif
        if i == 1
            gif(filename)
        else
            gif('DelayTime',1/60)
        end
    else
        pause(1/60)
    end
end

%% IMPORTING DATA FROM COMSOL

% import data from comsol
d=readmatrix('midpoint_final.txt','headerlines',9, 'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join');

plot(d(:,1), d(:,2))

% import data from comsol, again
d2=readmatrix('lineslice170.txt','headerlines',9, 'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join');
%% PLOTTING CENTRE POINT EVALUATION ANIMATED
figure()
size_file = size(d2);

for time = 1:size_file(2)
    if time == 1
        gif('comsolvmatlabanimatedleshgo.gif')
    else
        gif('DelayTime', 1/60)
    end
    plot(d2(:,1),d2(:,time))
    hold on
    midpoint = u(length(y)/2,:,time);
    plot(x, midpoint)
    legend('COMSOL','MATLAB')
    ylim([-0.8 1.5])
    title('Slice at y=y/2 for all x over time, MATLAB vs COMSOL')
    hold off    
end

%% PLOTTING CENTRE POINT EVALUATION
figure()
midpoints = zeros(1,n);
for time = 1:length(t)
    midpoints(time) = u(length(y)/2,length(x)/2,time);
end
plot(t,midpoints)
hold on
plot(d(:,1), d(:,2))
error = midpoints-d(:,2)
plot(t,error)
title('Point evaluation of MATLAB solution vs COMSOL solution')
xlabel('Time [s]')
ylabel('Z [m]')
legend('MATLAB','COMSOL','Error')


%% FUNCTION TO GENERATE GRID POINTS AND BIAS THEM TOWARDS EDGE OF DOMAIN
function grid_points = generate_grid(L,N,alpha)
% DO NOT
% PASS
% ODD NUMBERS
% INTO
% THIS FUNCTION
x_k = @(alpha,K,L,n) L*alpha*sum((1+alpha).^(K-1))/(2*((1+alpha)^n-1));

% CALCULATE OUT TO HALFWAY, THEN MIRROR AND COPY FOR THE OTHER HALF
halfway_point = N/2;
X = zeros(1,halfway_point);
for i  = 1:halfway_point-1
    K = 1:i;
    X(i+1) = x_k(alpha, K, L, halfway_point);
end

grid_points = zeros(1,N);
grid_points(1, 1:halfway_point) = X;
grid_points(1, halfway_point+1:N) = L-fliplr(X);
end


