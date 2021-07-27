function [K] = local_bar(E,A,L)
% returns stiffness vector K given E A and L

K = E*A/L * [1 -1; -1 1];
end

