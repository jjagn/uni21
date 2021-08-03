function [K_hat, lambda] = global_frame(K, alpha)
% Outputs stiffness matrix in global coords K hat and lamba of frame given 
% stiffness matrix K and angle to global coords alpha
c = cosd(alpha);
s = sind(alpha);

lambda33 = [c s 0;...
            -s c 0;...
            0 0 1];
        
zero_pad = zeros(3);

lambda = [lambda33 zero_pad;...
          zero_pad lambda33];
      
K_hat = lambda' * K * lambda;
end

