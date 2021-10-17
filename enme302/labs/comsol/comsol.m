close all; clc;

d=readmatrix('data1.txt','headerlines',9);
figure()
plot(d(:,1),d(:,3:2:end)-273.15);
figure()
plot(d(:,1),d(:,3:2:end)-d(:,4:2:end));