clear; close all; clc

fourier_terms = 5;

syms x y t
eqn1 = 0

for a = 1:fourier_terms
    for b = 1:fourier_terms
        eqn1 = eqn1 + ((1+(-1)^(a+1))*(1+(-1)^(b+1)))/(a^3*b^3)*sin(a*pi/2*x)*sin(b*pi/3*y)*cos(pi*sqrt(9*a^2+4*b^2)*t)
    end
end



