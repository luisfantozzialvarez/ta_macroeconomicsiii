syms gamma beta kappa delta,
A = [0 0 beta; -1/gamma 1 1/gamma; 1 0 0];
B = [0 -kappa 1; 0 1 0; 0 0 delta];
F = inv(A)*B;
e = eig(F)


A = [0 beta; 1 1/gamma];
B = [-kappa 1; 1 delta/gamma];
F = inv(A)*B;
e = eig(F)