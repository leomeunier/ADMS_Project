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

% Change the value of frequency range for the numerical identification 
Fmax = 200;
fs = 0.01:0.01:Fmax;

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
v = 0;

for k = gammasol
    v = v+1;
    Hcap = [1 0 1;  -sin(k*L) cosh(k*L) sinh(k*L) ; -cos(k*L) sinh(k*L) cosh(k*L)];    
    N = [0 ; -cos(k*L) ;sin(k*L)];
    Xicap = -inv(Hcap)*N;
    Xi = [1 Xicap.'];
    X(v,:)= Xi;
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

% 1. Getting the modal masses

damp = 1/100;
m1 = trapz(m*Phi1.^2)*damp;
m2 = trapz(m*Phi2.^2)*damp;
m3 = trapz(m*Phi3.^2)*damp;
m4 = trapz(m*Phi4.^2)*damp;

xk = 1.2;
xj = 0.2;

% 2. Getting the mode values in xj and xk

Phi1j = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xj);
Phi2j = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xj);
Phi3j = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xj);
Phi4j = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xj);
Phi1k = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xk);
Phi2k = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xk);
Phi3k = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xk);
Phi4k = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xk);

% 3. Calculating the FRF

Gjk1 = ((Phi1j*Phi1k/m1)./(-omega.^2 + 1i*2*damp*omegasol(1)*omega + omegasol(1)^2)) + ((Phi2j*Phi2k/m2)./(-omega.^2 + 1i*2*damp*omegasol(2)*omega + omegasol(2)^2)) + ...
    ((Phi3j*Phi3k/m3)./(-omega.^2 + 1i*2*damp*omegasol(3)*omega + omegasol(3)^2)) + ((Phi4j*Phi4k/m4)./(-omega.^2 + 1i*2*damp*omegasol(4)*omega + omegasol(4)^2));

% 4. Plot of the FRF (real then imaginary part)
figure(3)
subplot(2,1,1)
semilogy(fs,abs(Gjk1))
subplot(2,1,2)
plot(fs,angle(-Gjk1))


% Trying to estimate A1 for the least square minimization
estim1 = imag(Gjk1(4.5/0.01))*2*omegasol(1)^2*damp;

% Computing FRFs to get all experimental values
xk = 0.5;
xj = 1;

Phi1j = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xj);
Phi2j = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xj);
Phi3j = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xj);
Phi4j = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xj);
Phi1k = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xk);
Phi2k = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xk);
Phi3k = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xk);
Phi4k = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xk);

Gjk2 = ((Phi1j*Phi1k/m1)./(-omega.^2 + 1i*2*damp*omegasol(1)*omega + omegasol(1)^2)) + ((Phi2j*Phi2k/m2)./(-omega.^2 + 1i*2*damp*omegasol(2)*omega + omegasol(2)^2)) + ...
    ((Phi3j*Phi3k/m3)./(-omega.^2 + 1i*2*damp*omegasol(3)*omega + omegasol(3)^2)) + ((Phi4j*Phi4k/m4)./(-omega.^2 + 1i*2*damp*omegasol(4)*omega + omegasol(4)^2));

estim2 = (Phi1j*Phi1k/m1) + (Phi2j*Phi2k/m2) + (Phi3j*Phi3k/m3) + (Phi4j*Phi4k/m4);
xk = 0.1;
xj = 0.8;

Phi1j = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xj);
Phi2j = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xj);
Phi3j = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xj);
Phi4j = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xj);
Phi1k = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xk);
Phi2k = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xk);
Phi3k = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xk);
Phi4k = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xk);

Gjk3 = ((Phi1j*Phi1k/m1)./(-omega.^2 + 1i*2*damp*omegasol(1)*omega + omegasol(1)^2)) + ((Phi2j*Phi2k/m2)./(-omega.^2 + 1i*2*damp*omegasol(2)*omega + omegasol(2)^2)) + ...
    ((Phi3j*Phi3k/m3)./(-omega.^2 + 1i*2*damp*omegasol(3)*omega + omegasol(3)^2)) + ((Phi4j*Phi4k/m4)./(-omega.^2 + 1i*2*damp*omegasol(4)*omega + omegasol(4)^2));

estim3 = (Phi1j*Phi1k/m1) + (Phi2j*Phi2k/m2) + (Phi3j*Phi3k/m3) + (Phi4j*Phi4k/m4);
xk = 1;
xj = 0.4;

Phi1j = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xj);
Phi2j = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xj);
Phi3j = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xj);
Phi4j = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xj);
Phi1k = phi(X(1,1),X(1,2),X(1,3),X(1,4),gammasol(1),xk);
Phi2k = phi(X(2,1),X(2,2),X(2,3),X(2,4),gammasol(2),xk);
Phi3k = phi(X(3,1),X(3,2),X(3,3),X(3,4),gammasol(3),xk);
Phi4k = phi(X(4,1),X(4,2),X(4,3),X(4,4),gammasol(4),xk);

Gjk4 = ((Phi1j*Phi1k/m1)./(-omega.^2 + 1i*2*damp*omegasol(1)*omega + omegasol(1)^2)) + ((Phi2j*Phi2k/m2)./(-omega.^2 + 1i*2*damp*omegasol(2)*omega + omegasol(2)^2)) + ...
    ((Phi3j*Phi3k/m3)./(-omega.^2 + 1i*2*damp*omegasol(3)*omega + omegasol(3)^2)) + ((Phi4j*Phi4k/m4)./(-omega.^2 + 1i*2*damp*omegasol(4)*omega + omegasol(4)^2));

estim4 = (Phi1j*Phi1k/m1) + (Phi2j*Phi2k/m2) + (Phi3j*Phi3k/m3) + (Phi4j*Phi4k/m4);


% Plotting all the experimentals FRFs

figure(4)
semilogy(fs,abs(Gjk1))
hold on
semilogy(fs,abs(Gjk2))
hold on
semilogy(fs,abs(Gjk3))
hold on
semilogy(fs,abs(Gjk4))


% Building GrEXP matrix 

GrEXP = [ Gjk1 ; Gjk2 ; Gjk3 ; Gjk4 ];

% Least square minimization (have a look on epsilon function)



% Initial guesses
x0 = [0.1778,0.7369,0.0270,0.4946,1/100,1/100,1/100,1/100,4.5*2*pi,28.2*2*pi,79*pi*2,154.9*pi*2,0,0,0,0];
% Lsqm function
x = lsqnonlin(@(x) epsilon(x, fs, GrEXP), x0);
% Having a first look at x
disp(x);

% Computing of GrNUMs to compare them to Gjk

%  each GrNUMi is supposed to match with the peak number i 

GrNUM1 = (-x(4)./(-omega.^2 + 1i*2*x(8)*x(12)*omega + x(12)^2)) + x(13)./(omega.^2) + 1i*x(14)./(omega.^2) + x(15) + 1i*x(16);

figure(5)
subplot(2,1,1)
semilogy(fs,abs(Gjk1))
hold on 
semilogy(fs,abs(GrNUM1),'o')
subplot(2,1,2)
plot(fs,angle(Gjk1))
hold on 
plot(fs,angle(GrNUM1),'o')








