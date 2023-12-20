%% Initialisation of the problem

clear all
clc

[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure;

dis_stru(posit,l,gamma,xy,pr,idb,ndof); % for drawing the bike

omegamax = 200*2*pi;
omega1 = ((pi./l).^2).*sqrt(EJ./m);

Verif = omega1/omegamax;

N = nnod*3;

[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

case_studied = 0; % Use 0 if you study the inp file without A' (Q1,Q2,Q3,Q4 n Q7) , Use 1 for the inp file with A' (Q5 n Q6)

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

Eml1 = zeros(N,3);
Eml3 = zeros(N,3);
Eml8 = zeros(N,3);
Eml10 = zeros(N,3);

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

Ek1 = zeros(N,1);
Ek1(idb(1,2),1) = 1;
K1 = Ek1*k*(Ek1');
K = K + K1;

% this one is with ADMS_Assignement2.inp

if case_studied == 0
    Ek8 = zeros(45,1);
    Ek8(idb(8,2),1) = 1;
    K8 = Ek8*k*(Ek8');
    K = K + K8;
end

% this one is with ADMS_Assignement5.inp 

if case_studied==1
    K_k_L = [0 1 0 0 -1 0]'*k*[0 1 0 0 -1 0];
    g = 0;
    lambda_k = [cos(g) sin(g) 0
    -sin(g) cos(g) 0
    0 0 1];
    Lambda_k = [lambda_k zeros(3,3)
    zeros(3,3) lambda_k ];
    K_k_G = Lambda_k'* K_k_L * Lambda_k;
    idofn8 = idb(8,:);
    idofn16 = idb(16,:);
    idofk = [idofn8 idofn16];
    K(idofk, idofk) = K(idofk, idofk) + K_k_G;
end


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
Nforce = EA(n_el)*b;
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
semilogy(fs,abs(Nforce))
subplot(2,1,2)
plot(fs,angle(Nforce))

% Constraint force in C

ncons = nnod*3 - ndof;
for ii = 1:length(om)
R(:,ii)= (-MCF*om(ii)^2 + 1i*CCF*om(ii) + KCF)*X(:,ii);
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

figure(13)
subplot(2,1,1)
semilogy(fs,abs(-om.^2.*FRF_modvf))
hold on
semilogy(fs,abs(-om.^2.*X(vf,:)))
subplot(2,1,2)
plot(fs,angle(-om.^2.*FRF_modvf))
hold on
plot(fs,angle(-om.^2.*X(vf,:)))

figure(14)
subplot(2,1,1)
semilogy(fs,abs(-om.^2.*FRF_modhh))
hold on
semilogy(fs,abs(-om.^2.*X(hh,:)))
subplot(2,1,2)
plot(fs,angle(-om.^2.*FRF_modhh))
hold on
plot(fs,angle(-om.^2.*X(hh,:)))

% What we can notice in these comparisons is that the modal approach
% follows well the "experimental" curve for the two first peaks, that are
% related to the first two natural frequencies and so the first two modes. 

%% Input : vertical displacement of point A' Output : acceleration at node F (use ADMS_Assignement5.inp)

vf = idb(10,2);
fs = 0:0.1:200;
om = fs*2*pi;
ncons = nnod*3-ndof;
Y = zeros(ncons,1);
Y(idb(16,2)-ndof) = 1;
for ii = 1:length(om)
X(:,ii) = -(-MFF*om(ii)^2 + 1i*CFF*om(ii) + KFF)\ ...
(-MFC*om(ii)^2 + 1i*CFC*om(ii) + KFC)*Y; 
end

figure(15)
subplot(2,1,1)
semilogy(fs,abs(-(om.^2).*X(vf,:)))
subplot(2,1,2)
plot(fs,angle(-(om.^2).*X(vf,:)))

%% Time history of the steady-state vertical acceleration of point H 

% We divide the problem in 2 and then we will use the superposition
% approach ( Use also ADMS_Assignement5.inp here ) 

% For first irregularity

v = 12; % speed of bike
lamb1 = 1;
f1 = v/lamb1;
om1 = 2*pi*f1;
A1 = 0.001;
ncons = nnod*3-ndof;
Y1 = zeros(ncons,1);
Y1(idb(16,2)-ndof) = A1; %[m]
xx1 = -(-MFF*om1^2 + 1i*CFF*om1 + KFF)\ ...
(-MFC*om1^2 + 1i*CFC*om1 + KFC)*Y1; 
x1 = xx1(idb(13,2));

dt = 0.0001;
t = 0:dt:1;
x1_t = zeros(length(t),1);
a1_t = zeros(length(t),1);
for ii = 1:length(t)
x1_t(ii) = abs(x1)*cos(om1*t(ii)+angle(x1));
a1_t(ii) = -om1^2*abs(x1)*cos(om1*t(ii)+angle(x1));
end

% For second irregularity

lamb2 = 0.6;
f2 = v/lamb2;
om2 = 2*pi*f2;
A2 = 0.0005;
ncons = nnod*3-ndof;
Y2 = zeros(ncons,1);
Y2(idb(16,2)-ndof) = A2; %[m]
xx2 = -(-MFF*om2^2 + 1i*CFF*om2 + KFF)\ ...
(-MFC*om2^2 + 1i*CFC*om2 + KFC)*Y2; 
x2 = xx2(idb(13,2));

dt = 0.0001;
t = 0:dt:1;
x2_t = zeros(length(t),1);
a2_t = zeros(length(t),1);
for ii = 1:length(t)
x2_t(ii) = abs(x2)*cos(om2*t(ii)+angle(x2));
a2_t(ii) = -om2^2*abs(x2)*cos(om2*t(ii)+angle(x2));
end

% Superposition 

a_t = a1_t + a2_t;
figure(16)
plot(t,a_t)

%% Static response of the structure due to the weight of the cyclist

FG  = zeros(N,1);
indexH = idb(13,2);
indexF = idb(10,2);
FG(indexH) = -600;
FG(indexF) = -100;
xF = KFF \ FG(1:ndof);
figure(17)
diseg2(xF,100,incid,l,gamma,posit,idb,xy)
title(['static deflection']);

disp(xF(idb(11,2)))


