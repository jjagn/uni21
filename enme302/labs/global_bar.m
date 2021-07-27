function [Khat,lambda] = global_bar(K,alpha)
% takes stiffness matrix K and angle alpha and outputs K hat and lambda 
% matrices
c = cosd(alpha);
s = sind(alpha);

lambda = [c s 0 0; 0 0 c s];

Khat = lambda' * K * lambda;

end

