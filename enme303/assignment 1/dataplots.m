clear
clc
close all

p1 = import_enme303('jdw144_p1');
pd1 = import_enme303('jdw144_pd1');
pd2 = import_enme303('jdw144_pd2');
pid1 = import_enme303('jdw144_pid1');

figure()
plot(p1(:, 1), p1(:, [2,4]));
title('Cart position and commanded position against time for P controller 1')
xlim([0, 6])
legend('cart position', 'position command')

figure()
plot(pd1(:, 1), pd1(:, [2,4]));
title('Cart position and commanded position against time for PD controller 1')
xlim([0, 6])
legend('cart position', 'position command')

figure()
plot(pd2(:, 1), pd2(:, [2,4]));
title('Cart position and commanded position against time for PD controller 2')
xlim([0, 5])
legend('cart position', 'position command')

figure()
plot(pid1(:, 1), pid1(:, [2,4]));
title('Cart position and commanded position against time for PID controller 1')
xlim([0, 17.9])
legend('cart position', 'position command')