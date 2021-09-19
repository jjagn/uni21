clear; close all; clc

% [tip_x, tip_y, tip_rot, upper_x, upper_y, upper_z, lower_x, lower_y,
% lower_z, maxmimum_total_normal_stress]

original = [83.33 25 0 83.33 0 0];
new = [83.08 24.8 0.0619 83.08 0.201 0.2397];


delta = original - new
delta_pct = delta./original .* 100

