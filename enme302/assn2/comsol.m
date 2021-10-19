clear; close all; clc;

d=readmatrix('midpoint.txt','headerlines',9, 'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join');
d = d(3:end)
