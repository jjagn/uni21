function plot_deflected_shape(node1XG, node1YG, node2XG, node2YG, d_e, magFactor, n, L, alpha)
% plots deflected shape of beam using shape functions
x_e = linspace(0, L, n);
phi_1 = 1 - x_e ./ L;
phi_2 = x_e ./ L;

N1 = 1 - 3 .* x_e .^ 2 ./ L^2 + 2 .* x_e .^ 3 ./L .^3;
N2 = x_e .^ 3 ./ L .^ 2 - 2 .* x_e .^ 2 ./ L + x_e;
N3 = 3 .* x_e .^ 2 ./ L .^ 2 - 2 .* x_e .^3 ./ L .^3;
N4 = x_e .^ 3 ./ L .^ 2 - x_e .^ 2 ./ L;

u = phi_1 .* d_e(1) + phi_2 .* d_e(4);
v = N1 * d_e(2) +  N2 * d_e(3) +  N3  * d_e(5) + N4 * d_e(6);

deflections_XG = u .* cosd(alpha) - v .* sind(alpha);
deflections_YG = u .* sind(alpha) + v .* cosd(alpha);

undeflected_XG = linspace(node1XG, node2XG, n);
undeflected_YG = linspace(node1YG, node2YG, n);

deflected_XG = undeflected_XG + magFactor .* deflections_XG;
deflected_YG = undeflected_YG + magFactor .* deflections_YG;

plot(undeflected_XG, undeflected_YG, 'b-')
hold on
plot(deflected_XG, deflected_YG, 'r-')
legend('original shape', 'deflected shape')
end

