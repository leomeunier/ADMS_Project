clear all 

% Variables

L = 1200e-3;
h = 8e-3;
b = 40e-3;
rho = 2700;
E = 68e9;
J = (1/3)*b*h^3;
m = rho*b*h*L;
p = pi;

% Identification des frequences propres

Fmax = 200;
fs = 0:0.05:Fmax;

omega = 2*pi*fs;

gamma = (m/(E*J)^(1/4))*sqrt(omega);

det = (2 + 2*cos(gamma*L).*cosh(gamma*L)); %.*(gamma.^6);

%semilogy(fs, abs(det))

% Identification des modes propres

fsol = [4.5 28 69 155];
omegasol = 2*p*fsol;
gammasol = (m/(E*J)^(1/4))*sqrt(omegasol);
X = zeros(4,4);
i = 0;

for k = gammasol;
    i = i+1;
    Hcap = [1 0 1;  -sin(k*L) cosh(k*L) sinh(k*L) ; -cos(k*L) sinh(k*L) cosh(k*L)];    
    N = [0 ; -cos(k*L) ;sin(k*L)];
    Xicap = -inv(Hcap)*N;
    Xi = [1 Xicap.'];
    X(i,:)= Xi;
end

x = 0:0.1:1.2;

Phi1 = X(1,1)*sin(gammasol(1)*x) + X(1,2)*cos(gammasol(1)*x) + X(1,3)*sinh(gammasol(1)*x) + X(1,4)*cosh(gammasol(1)*x);
Phi2 = X(2,1)*sin(gammasol(2)*x) + X(2,2)*cos(gammasol(2)*x) + X(2,3)*sinh(gammasol(2)*x) + X(2,4)*cosh(gammasol(1)*x);
Phi3 = X(3,1)*sin(gammasol(3)*x) + X(3,2)*cos(gammasol(3)*x) + X(3,3)*sinh(gammasol(3)*x) + X(3,4)*cosh(gammasol(1)*x);
Phi4 = X(4,1)*sin(gammasol(4)*x) + X(4,2)*cos(gammasol(4)*x) + X(4,3)*sinh(gammasol(4)*x) + X(4,4)*cosh(gammasol(1)*x);

%plot(x,Phi1./max(abs(Phi1)))
%plot(x,Phi2./max(Phi2))
%plot(x,Phi3./max(Phi3))
%plot(x,Phi4./max(Phi4))

% Computing the FRF

damp = 1/100;
xk = x(3);
xj = x(8);

m1 = trapz(m*Phi1.^2);
m2 = trapz(m*Phi2.^2);
m3 = trapz(m*Phi3.^2);
m4 = trapz(m*Phi4.^2);

Phi1j = X(1,1)*sin(gammasol(1)*xj) + X(1,2)*cos(gammasol(1)*xj) + X(1,3)*sinh(gammasol(1)*xj) + X(1,4)*cosh(gammasol(1)*xj);
Phi2j = X(2,1)*sin(gammasol(2)*xj) + X(2,2)*cos(gammasol(2)*xj) + X(2,3)*sinh(gammasol(2)*xj) + X(2,4)*cosh(gammasol(1)*xj);
Phi3j = X(3,1)*sin(gammasol(3)*xj) + X(3,2)*cos(gammasol(3)*xj) + X(3,3)*sinh(gammasol(3)*xj) + X(3,4)*cosh(gammasol(1)*xj);
Phi4j = X(4,1)*sin(gammasol(4)*xj) + X(4,2)*cos(gammasol(4)*xj) + X(4,3)*sinh(gammasol(4)*xj) + X(4,4)*cosh(gammasol(1)*xj);
Phi1k = X(1,1)*sin(gammasol(1)*xk) + X(1,2)*cos(gammasol(1)*xk) + X(1,3)*sinh(gammasol(1)*xk) + X(1,4)*cosh(gammasol(1)*xk);
Phi2k = X(2,1)*sin(gammasol(2)*xk) + X(2,2)*cos(gammasol(2)*xk) + X(2,3)*sinh(gammasol(2)*xk) + X(2,4)*cosh(gammasol(1)*xk);
Phi3k = X(3,1)*sin(gammasol(3)*xk) + X(3,2)*cos(gammasol(3)*xk) + X(3,3)*sinh(gammasol(3)*xk) + X(3,4)*cosh(gammasol(1)*xk);
Phi4k = X(4,1)*sin(gammasol(4)*xk) + X(4,2)*cos(gammasol(4)*xk) + X(4,3)*sinh(gammasol(4)*xk) + X(4,4)*cosh(gammasol(1)*xk);

Gjk = ((Phi1j*Phi1k/m1)./(-fs.^2 + i*2*damp*gammasol(1)*fs + gammasol(1)^2)) + ((Phi2j*Phi2k/m1)./(-fs.^2 + i*2*damp*gammasol(2)*fs + gammasol(2)^2)) + ...
    ((Phi3j*Phi3k/m1)./(-fs.^2 + i*2*damp*gammasol(3)*fs + gammasol(3)^2)) + ((Phi4j*Phi4k/m1)./(-fs.^2 + i*2*damp*gammasol(4)*fs + gammasol(4)^2));

%plot(fs,abs(Gjk))
%plot(fs,angle(Gjk))







