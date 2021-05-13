clear
clc
close all

p1 = import_enme303('jdw144_p1');
pd1 = import_enme303('jdw144_pd1');
pd2 = import_enme303('jdw144_pd2');
pid1 = import_enme303('jdw144_pid1');

figure('DefaultAxesFontSize', 14)
plot(p1(:, 1), p1(:, [2,4]), 'LineWidth', 1);
grid ON
xlim([0, 6])
legend('Cart position', 'Position command')
ylabel('Position (m)')
xlabel('Time (s)')

figure('DefaultAxesFontSize', 14)
plot(pd1(:, 1), pd1(:, [2,4]), 'LineWidth', 1);
grid ON
xlim([0, 10])
ylim([-0.11, 0.11])
ylabel('Position (m)')
xlabel('Time (s)')
legend('Cart position', 'Position command')

% figure()
% plot(pd2(:, 1), pd2(:, [2,4]), 'LineWidth', 1);
% grid ON
% title('Cart position and commanded position against time for PD controller 2')
% xlim([0, 5])
% legend('cart position', 'position command')

figure('DefaultAxesFontSize', 14)
plot(pid1(:, 1), pid1(:, [2,4]), 'LineWidth', 1);
grid ON
xlim([0, 17.9])
legend('Cart position', 'Position command')
ylabel('Position (m)')
xlabel('Time (s)')