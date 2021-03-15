% num = [3 9];
% den = [1 2 9];
% 
% figure();
% impulse(num, den);
% figure();
% step(num, den);
% figure();
% pzmap(num, den);
% 
% numH = 160*conv([1 2 5], [1 0 7]);
% denH = conv([1 5 40], [1 0.03 0.06]);
% 
% figure();
% damp(denH)
% 
% num1 = 1;
% den1 = conv([1 0], conv([1 1], [1 2]));
% 
rlocus(num1, den1);
title('root locus');

[K, POLES] = rlocfind(num1, den1);

figure();
num = K;
den3 = conv([1 0], conv([1 1], [1 2]));
[numcl, dencl] = feedback(num, den, 1, 1);
impulse(numcl, dencl);
title('closed loop system response');