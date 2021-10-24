clear; close all; clc
alpha = 0.05;
L = 2;
H = 3;
N = 10;

x_points = generate_grid(L,N*L,alpha);
y_points = generate_grid(H,N*H,alpha);

[X,Y] = meshgrid(x_points, y_points);

surf(X,Y,X*0)

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