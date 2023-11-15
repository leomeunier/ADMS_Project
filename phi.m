function [Phi] = phi(A,B,C,D,gamma,x)

Phi = A*sin(gamma*x) + B*cos(gamma*x) + C*sinh(gamma*x) + D*cosh(gamma*x);

end
