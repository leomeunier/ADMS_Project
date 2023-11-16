function [Phi] = phi(A,B,C,D,gamma,x)

Phi = A*cos(gamma*x) + B*sin(gamma*x) + C*cosh(gamma*x) + D*sinh(gamma*x);

end
