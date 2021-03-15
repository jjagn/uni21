num = 1;
den = [1 0.2 1];
bode(num, den);
title('Bode plot');
figure()

axis([-2 2 -2 2]);
nyquist(num, den);