function plotDeflectedBar(node1XG, node1YG, D_e, magFactor, n, L, alpha)
% plots deflected shape of beam using shape functions

node2XG = node1XG + L * cosd(alpha);
node2YG = node1YG + L * sind(alpha);

node1DeflectedXG = node1XG + D_e(1) * magFactor;
node1DeflectedYG = node1YG + D_e(2) * magFactor;
node2DeflectedXG = node2XG + D_e(3) * magFactor;
node2DeflectedYG = node2YG + D_e(4) * magFactor;

undeflected_XG = linspace(node1XG, node2XG, n);
undeflected_YG = linspace(node1YG, node2YG, n);

deflected_XG = linspace(node1DeflectedXG, node2DeflectedXG, n);
deflected_YG = linspace(node1DeflectedYG, node2DeflectedYG, n);

plot(undeflected_XG, undeflected_YG, 'b-')
hold on
plot(deflected_XG, deflected_YG, 'r-')
legend('original shape', 'deflected shape')
end

