clear all 

% Variables

L = 1200e-3;
h = 8e-3;
b = 40e-3;
rho = 2700;
E = 68e9;
J = (1/12)*b*h^3;
m = rho*b*h;
p = pi;

% Identification of natural frequencies

Fmax = 200;
fs = 0:0.5:Fmax;

omega = 2*pi*fs;

gamma = ((m*omega.^2)/(E*J)).^(1/4);

det = (2 + 2*cos(gamma*L).*cosh(gamma*L)).*(gamma.^6);

figure(1)
semilogy(fs, abs(det))

% Identification natural modes

fsol = [4.5 28.2 79 154.9];
omegasol = 2*p*fsol;
gammasol = ((m*omegasol.^2)/(E*J)).^(1/4);
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

x = 0:0.01:1.2;

Phi1 = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),x);
Phi2 = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),x);
Phi3 = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),x);
Phi4 = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),x);


figure(2)
subplot(4,1,1)
plot(x,Phi1./max(abs(Phi1)))
subplot(4,1,2)
plot(x,Phi2./max(abs(Phi2)))
subplot(4,1,3)
plot(x,Phi3./max(abs(Phi3)))
subplot(4,1,4)
plot(x,Phi4./max(abs(Phi4)))

% Computing the FRF

damp = 1/100;
xk = 0.2;
xj = 1.2;

m1 = trapz(m*Phi1.^2);
m2 = trapz(m*Phi2.^2);
m3 = trapz(m*Phi3.^2);
m4 = trapz(m*Phi4.^2);

Phi1j = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xj);
Phi2j = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xj);
Phi3j = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xj);
Phi4j = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xj);
Phi1k = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xk);
Phi2k = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xk);
Phi3k = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xk);
Phi4k = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xk);

Gjk = ((Phi1j*Phi1k/m1)./(-omega.^2 + i*2*damp*omegasol(1)*omega + omegasol(1)^2)) + ((Phi2j*Phi2k/m2)./(-omega.^2 + i*2*damp*omegasol(2)*omega + omegasol(2)^2)) + ...
    ((Phi3j*Phi3k/m3)./(-omega.^2 + i*2*damp*omegasol(3)*fs + omegasol(3)^2)) + ((Phi4j*Phi4k/m4)./(-omega.^2 + i*2*damp*omegasol(4)*fs + omegasol(4)^2));

figure(3)
subplot(2,1,1)
semilogy(fs,abs(Gjk))
subplot(2,1,2)
plot(fs,angle(Gjk))









