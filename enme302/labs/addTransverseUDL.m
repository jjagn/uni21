function [f_eq] = addVerticalUDL(w_bar, L)
% creates f vector for adding UDL to frame element from load intensity w
% bar
f_eq = [0;...
        w_bar * L/2;...
        w_bar * L^2/12;...
        0;...
        w_bar * L/2;...
        w_bar * L^2/12];
end
