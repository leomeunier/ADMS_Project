%% Initialisation of the problem

clear all
clc

[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure;

dis_stru(posit,l,gamma,xy,pr,idb,ndof); % for drawing the bike

omegamax = 200*2*pi;
omega1 = ((pi./l).^2).*sqrt(EJ./m);

Verif = omega1/omegamax;

[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

%% Added Masses

% Masses in each node

m1 = 2;
m3 = 2;
m8 = 1;
m10 = 1.5;

% Matrix linked to these masses (No J because no radius)

M1c = [m1 0 0; 0 m1 0; 0 0 0];
M3c = [m3 0 0; 0 m3 0; 0 0 0];
M8c = [m8 0 0; 0 m8 0; 0 0 0];
M10c = [m10 0 0; 0 m10 0; 0 0 0];

% Construction of the full mass matrices

Eml1 = zeros(45,3);
Eml3 = zeros(45,3);
Eml8 = zeros(45,3);
Eml10 = zeros(45,3);

Eml1(idb(1,1),1)= 1;
Eml1(idb(1,2),2)= 1;
Eml1(idb(1,3),3)= 1;

Eml3(idb(3,1),1)= 1;
Eml3(idb(3,2),2)= 1;
Eml3(idb(3,2),3)= 1;

Eml8(idb(8,1),1)= 1;
Eml8(idb(8,1),2)= 1;
Eml8(idb(8,1),3)= 1;

Eml10(idb(10,1),1)= 1;
Eml10(idb(10,2),2)= 1;
Eml10(idb(10,3),3)= 1;

% Full mass matrix

M1 = Eml1*M1c*(Eml1');
M3 = Eml3*M3c*(Eml3');
M8 = Eml8*M8c*(Eml8');
M10 = Eml10*M10c*(Eml10');


M = M + M1 + M3 + M8 + M10;

%% Added Springs 

k = 4e+05;

Ek1 = zeros(45,1);
Ek8 = zeros(45,1);

Ek1(idb(1,2),1) = 1;
Ek8(idb(8,2),1) = 1;

K1 = Ek1*k*(Ek1');
K8 = Ek8*k*(Ek8');

% final spring matrix

K = K + K1 + K8;


%% Getting the natural frequencies and mode shapes

MFF = M(1:ndof,1:ndof);
KFF = K(1:ndof,1:ndof);
MFC = M(1:ndof,ndof+1:nnod*3);
KFC = K(1:ndof,ndof+1:nnod*3);
MCF = M(ndof+1:nnod*3,1:ndof);
KCF = K(ndof+1:nnod*3,1:ndof);
MCC = M(ndof+1:nnod*3,ndof+1:nnod*3);
KCC = K(ndof+1:nnod*3,ndof+1:nnod*3);

[modes, omega2] = eig(inv(MFF)*KFF);
omega = diag(sqrt(omega2));
% Sort frequencies in ascending order
[omega_sorted, omega_sorted_indices] = sort(omega);
% Sort mode shapes in ascending order
modes_sorted = modes(:,omega_sorted_indices);

disp(omega_sorted/(2*pi));

%% Drawing modes shapes

figure(2);
scale_factor = 2;
for ii = 1:4
    mode = modes_sorted(:,ii);
    subplot(2,2,ii)
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy);
end

%% Damping matrix 

alpha = 6;
beta = 10e-5;
CFF = alpha*MFF + beta*KFF;
CFC = alpha*MFC + beta*KFC;
CCF = alpha*MCF + beta*KCF;
CCC = alpha*MCC + beta*KCC;


%% Frequency Response Function

% Get the value in idb of C ( 3rd node so 3rd row ) for Y displacement (
% 2nd column)

F0 = zeros(ndof,1);
index = idb(3,2);
F0(index) = 1;
fs = 0:0.1:200;
om = fs*2*pi;
for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

% Vertical displacement of node F

vf = idb(10,2);

figure(3)
subplot(2,1,1)
semilogy(fs,abs(X(vf,:)))
subplot(2,1,2)
plot(fs,angle(X(vf,:)))

% Vertical acceleration of node F 

figure(4)
subplot(2,1,1)
semilogy(fs,abs(-(om.^2).*X(vf,:)))
subplot(2,1,2)
plot(fs,angle(-(om.^2).*X(vf,:)))

%  Horizontal displacement of node H

hh = idb(13,1);

figure(5)
subplot(2,1,1)
semilogy(fs,abs(X(hh,:)))
subplot(2,1,2)
plot(fs,angle(X(hh,:)))

% Horizontal acceleration of node H 

figure(6)
subplot(2,1,1)
semilogy(fs,abs(-(om.^2).*X(hh,:)))
subplot(2,1,2)
plot(fs,angle(-(om.^2).*X(hh,:)))


% Shear force T, bending moment M and axial force N evaluated in the midpoint of the GE tube 

mge = idb(11,:);

n_el = 11; 
L_el = l(n_el);
idof_i = idb(11,:); 
idof_j = idb(12,:); 
lambda = [cos(gamma(n_el)) sin(gamma(n_el)) 0; 
-sin(gamma(n_el)) cos(gamma(n_el)) 0; 
0 0 1]; 
Xi = lambda*X(idof_i,:);
Xj = lambda*X(idof_j,:); 
b = Xi(3,:);
c = -3/L_el^2*Xi(2,:) +3/L_el^2*Xj(2,:) -2/L_el^1*Xi(3,:) -1/L_el^1*Xj(3,:);
d = 2/L_el^3*Xi(2,:) -2/L_el^3*Xj(2,:) +1/L_el^2*Xi(3,:) +1/L_el^2*Xj(3,:);
T = EJ(n_el)*(6*d);
N = EA(n_el)*b;
Mbend = EJ(n_el)*(2*c + 6*d*0);

figure(7)
subplot(2,1,1)
semilogy(fs,abs(T))
subplot(2,1,2)
plot(fs,angle(T))

figure(8)
subplot(2,1,1)
semilogy(fs,abs(Mbend))
subplot(2,1,2)
plot(fs,angle(Mbend))

figure(9)
subplot(2,1,1)
semilogy(fs,abs(N))
subplot(2,1,2)
plot(fs,angle(N))

% Constraint force in C

ncons = nnod*3 - ndof;
Y = zeros(ncons,1);
Y(idb(3,1)-ndof) = 0.01;
for ii = 1:length(om)
Xc(:,ii)= -(-MFF*om(ii)^2 + 1i*CFF*om(ii) + KFF)\...
(-MFC*om(ii)^2 + 1i*CFC*om(ii) + KFC)*Y; 
R(:,ii)= (-MCF*om(ii)^2 + 1i*CCF*om(ii) + KCF)*Xc(:,ii)+...
(-MCC*om(ii)^2 + 1i*CCC*om(ii) + KCC)*Y; 
end
RxC = R(idb(3,1)-ndof,:);

figure(10)
subplot(2,1,1)
semilogy(fs,abs(RxC))
subplot(2,1,2)
plot(fs,angle(RxC))

%% Modal superposition approach

% Modal matrices
ii = 1:2; % first 2 mode shapes
Phi = modes_sorted(:,ii); 
Mmod = Phi'*MFF*Phi; 
Kmod = Phi'*KFF*Phi; 
Cmod = Phi'*CFF*Phi; 
Fmod = Phi'*F0;
% FRF in modal superposition approach
for ii = 1:length(om)
xx_mod(:,ii) = (-om(ii)^2*Mmod + 1i*om(ii)*Cmod + Kmod) \ Fmod;
end
xx_2m = Phi * xx_mod; 

FRF_modvf = xx_2m(vf,:); % FRF of the vertical displacement of F 
FRF_modhh = xx_2m(hh,:); % FRF of the horizontal displacement of H

figure(11)
subplot(2,1,1)
semilogy(fs,abs(FRF_modvf))
hold on
semilogy(fs,abs(X(vf,:)))
subplot(2,1,2)
plot(fs,angle(FRF_modvf))
hold on
plot(fs,angle(X(vf,:)))

figure(12)
subplot(2,1,1)
semilogy(fs,abs(FRF_modhh))
hold on
semilogy(fs,abs(X(hh,:)))
subplot(2,1,2)
plot(fs,angle(FRF_modhh))
hold on
plot(fs,angle(X(hh,:)))
