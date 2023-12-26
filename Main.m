%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical System Dynamics
% FEM script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load input file and assemble structure
[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,T,posit,nbeam,pr]=loadstructure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble mass and stiffness matricies
[M,K] = assem(incid,l,m,EA,EJ,T,gamma,idb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add concentrated elements

% Concentrated spring in node 1
k = 4e+05;
i_ndof_spring = idb(1,2);
K(i_ndof_spring,i_ndof_spring) = K(i_ndof_spring,i_ndof_spring) + k; 

% Concentrated spring in node 2
k = 4e+05;
i_ndof_spring = idb(2,2);
K(i_ndof_spring,i_ndof_spring) = K(i_ndof_spring,i_ndof_spring) + k; 

% Concentrated masses
m1 = 1; %[kg]
m2 = 2; %[kg]
m3 = 2; %[kg]
m4 = 1.5; %[kg]

% m1
M1_h = diag([m1 m1 0]);
i_ndof_m1 = idb(1,:);
M(i_ndof_m1,i_ndof_m1) = M(i_ndof_m1,i_ndof_m1) + M1_h;

% m2
M2_h = diag([m2 m2 0]);
i_ndof_m2 = idb(3,:);
M(i_ndof_m2,i_ndof_m2) = M(i_ndof_m2,i_ndof_m2) + M2_h;

% m3
M3_h = diag([m3 m3 0]);
i_ndof_m3 = idb(2,:);
M(i_ndof_m3,i_ndof_m3) = M(i_ndof_m3,i_ndof_m3) + M3_h;

% m4
M4_h = diag([m4 m4 0]);
i_ndof_m4 = idb(6,:);
M(i_ndof_m4,i_ndof_m4) = M(i_ndof_m4,i_ndof_m4) + M4_h;

%% Compute natural frequencies and mode shapes
MFF = M(1:ndof,1:ndof);
MCF = M(ndof+1:end,1:ndof);
MFC = M(1:ndof,ndof+1:end);
MCC = M(ndof+1:end,ndof+1:end);

KFF = K(1:ndof,1:ndof);
KCF = K(ndof+1:end,1:ndof);
KFC = K(1:ndof,ndof+1:end);
KCC = K(ndof+1:end,ndof+1:end);

[modes, omega] = eig(inv(MFF)*KFF);
omega = sqrt(diag(omega));
[omega,i_omega] = sort(omega);
freq0 = omega/2/pi;
modes = modes(:,i_omega);

figure(2)
nmodes = 4;
scale_factor = 2;
for ii = 1:nmodes
    mode = modes(:,ii);
    subplot(2,2,ii)
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy)
    legend('Undeformed structure',['Mode ',num2str(ii),' - freq ',num2str(freq0(ii)),' Hz']);
end

%% Damping Matrix

alpha = 6;      % [s^-1]
beta = 1e-5;    % [s]

C = alpha*M + beta*K;
CFF = C(1:ndof, 1:ndof);
CCF = C(ndof+1:end, 1:ndof);
CFC = C(1:ndof, ndof+1:end);

%% 3) Frequency Response Function Fc

%% 3A) FRF of displacement and acceleration

% Fc related to vertical displacement and acceleration of node F

freq = 0:0.1:200;
Om = 2*pi*freq;
F0 = zeros(ndof,1);

idof_j = idb(6,2); % vertical displacement of node F
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = (-Om(ii)^2)*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end
FRF_nFv = X(idof_j,:);

figure(3)
subplot(2,1,1)
semilogy(freq, abs(FRF_nFv),'LineWidth',1.5)
grid on
xlabel('frequency [Hz]')
ylabel('|FRF|')
title('Vertical displacement of node F due to the force applied in node C')
subplot(2,1,2)
plot(freq, angle(FRF_nFv),'LineWidth',1.5)
grid on
xlabel('frequency [Hz]')
ylabel('\phi(FRF)')

% Acceleration of node F
for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
    acc(:,ii) = -Om(ii)^2*X(:,ii);
end
FRF_nFa = acc(idof_j,:);

figure(7)
subplot(2,1,1)
semilogy(freq, abs(FRF_nFa),'LineWidth',1.5);
grid on
xlabel('frequency [Hz]')
ylabel('|FRF|')
title('Vertical acceleration of node F due to the force applied in node C')
subplot(2,1,2)
plot(freq, angle(FRF_nFa),'LineWidth',1.5);
grid on
xlabel('frequency [Hz]')
ylabel('\phi(FRF)')


% Horizontal displacement and acceleration of node H

F0 = zeros(ndof,1);

idof_j = idb(8,1); % horizontal displacement of node H
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end
FRF_nHh = X(idof_j,:);

figure(8)
hold on 
subplot(2,1,1);
semilogy(freq, abs(FRF_nHh),'LineWidth',1.5);
grid on
xlabel('frequency [Hz]')
ylabel('|FRF|')
title('Horizontal displacement of node H due to the force applied in node C')
subplot(2,1,2);
plot(freq, angle(FRF_nHh),'LineWidth',1.5);
grid on
xlabel('frequency [Hz]')
ylabel('\phi(FRF)')

% Acceleration of node H (horizontal displacement of node H is known)
acc(:,ii) = -Om(ii)^2*X(:,ii);
FRF_nHa = acc(idof_j,:);

figure(9)
subplot(2,1,1)
semilogy(freq, abs(FRF_nHa), 'LineWidth',1.5);
grid on
xlabel('frequency [Hz]')
ylabel('|FRF|') 
title('Horizontal acceleration of node H due to the force applied in node C')
subplot(2,1,2)
plot(freq, angle(FRF_nHa),'LineWidth',1.5);
grid on
xlabel('frequency [Hz]')
ylabel('\phi(FRF)')

%% 3B) Internal forces in node 11
nel = 15;
F0 = zeros(ndof,1);

idof_12 = idb(12,:); 
idof_11 = idb(11,:);
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end
X_12 = X(idof_12,:);
X_11 = X(idof_11,:);

lambda = [cos(gamma(nel)) sin(gamma(nel)) 0;
        -sin(gamma(nel)) cos(gamma(nel)) 0;
         0                0               1];

L_el = l(nel);
csi = L_el;
loc_X_12 = lambda*X_12; %Xi
loc_X_11 = lambda*X_11; %Xj

% b for axial force N
b = (loc_X_11(1,:) - loc_X_12(1,:))/L_el;

N_11 = EA(nel)*b;

% c and d for shear force T and bending moment M
c = -3/L_el^2*loc_X_12(2,:) + 3/L_el^2*loc_X_11(2,:) - 2/L_el*loc_X_12(3,:) - 1/L_el*loc_X_11(3,:);
d = 2/L_el^3*loc_X_12(2,:) - 2/L_el^3*loc_X_11(2,:) + 1/L_el^2*loc_X_12(3,:) + 1/L_el^2*loc_X_11(3,:);

T_11 = EJ(nel)*(6*d);
M_11 = EJ(nel)*(2*c + 6*d*csi);

figure(10)
subplot(2,1,1)
semilogy(freq,abs(N_11), 'LineWidth',1.5)
hold on
grid on
xlabel('frequency [Hz]')
ylabel('Magnitude')
title(['Magnitude of the internal forces at mid point of GE tube'])
semilogy(freq,abs(T_11), 'LineWidth',1.5)
semilogy(freq,abs(M_11), 'LineWidth',1.5)
legend('Axial force', 'Shear force', 'Bending moment')

subplot(2,1,2)
plot(freq,angle(N_11), 'LineWidth',1.5)
hold on
grid on
xlabel('frequency [Hz]')
ylabel('Phase [rad]')
title(['Phase of the internal forces at mid point of GE tube'])
plot(freq,angle(T_11), 'LineWidth',1.5)
plot(freq,angle(M_11), 'LineWidth',1.5)
legend('Axial force', 'Shear force', 'Bending moment')


%% 3C) Reaction forces on the constraint

F0 = zeros(ndof,1);
idof_k = idb(3,2); % force in node C
F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    Ar = -Om(ii)^2*MCF + 1i*Om(ii)*CCF + KCF;
    X(:,ii) = A\F0;
    rr(:,ii) = Ar*X(:,ii);
end

ndof_RV = idb(3,1) - ndof;
FRF_RV  = rr(ndof_RV,:);

figure(11)
subplot(2,1,1)
semilogy(freq, abs(FRF_RV),'LineWidth',1.5);
xlabel('frequency [Hz]')
ylabel('Magnitude')
title('Magnitude of the constraint force of node C')
grid on
subplot(2,1,2)
plot(freq, angle(FRF_RV),'LineWidth',1.5);
grid on
xlabel('frequency [Hz]')
ylabel('Phase [rad]')
title('Phase of the constraint force of node C')

%% 4) Modal superposition approach

F0 = zeros(ndof,1);
idof_k = idb(3,2); % force in node C
F0(idof_k) = 1;

modes_2 = modes(:,1:2);

Mmod = modes_2'*MFF*modes_2;
Kmod = modes_2'*KFF*modes_2;
Cmod = modes_2'*CFF*modes_2;
Fmod = modes_2'*F0;

for ii = 1:length(Om)
    xx_mod(:,ii) = (-Om(ii)^2*Mmod + 1i*Om(ii)*Cmod + Kmod)\Fmod;
    acc_mod(:,ii) = -Om(ii)^2*xx_mod(:,ii);
end
xx_m = modes_2*xx_mod;
acc_m = modes_2*acc_mod;
FRFm_nFv = xx_m(idb(6,2),:);
FRFm_nFa = acc_m(idb(6,2),:);
FRFm_nHh = xx_m(idb(8,1),:);
FRFm_nHa = acc_m(idb(8,1),:);

% NODE F

figure(12)
subplot(2,1,1)
semilogy(freq,abs(FRF_nFv),'LineWidth',1.5)
grid on
hold on
semilogy(freq,abs(FRF_nFa),'LineWidth',1.5)
hold on
semilogy(freq,abs(FRFm_nFv),'LineWidth',1.5)
hold on
semilogy(freq,abs(FRFm_nFa),'LineWidth',1.5)
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Magnitude of the vertical displacement and acceleration of node F')
legend('Displacement of node F using FE method', 'Acceleration of node F using FE method', 'Displacement of node F using modal approach (first 2 natural freq.)', 'Acceleration of node F using modal approach (first 2 natural freq.)')

subplot(2,1,2)
plot(freq,angle(FRF_nFv),'LineWidth',1.5)
hold on
plot(freq,angle(FRF_nFa),'LineWidth',1.5)
grid on
hold on
plot(freq,angle(FRFm_nFv),'LineWidth',1.5)
hold on
plot(freq,angle(FRFm_nFa),'LineWidth',1.5)
xlabel('Frequency [Hz]');
ylabel('Phase[rad]');
title('Phase of the vertical displacement and acceleration of node F')
legend('Displacement of node F using FE method', 'Acceleration of node F using FE method', 'Displacement of node F using modal approach (first 2 natural freq.)', 'Acceleration of node F using modal approach (first 2 natural freq.)')


% NODE H
figure(13)
subplot(2,1,1)
grid on
semilogy(freq,abs(FRF_nHh),'LineWidth',1.5)
hold on
grid on
semilogy(freq,abs(FRF_nHa),'LineWidth',1.5)
hold on
semilogy(freq,abs(FRFm_nHh),'LineWidth',1.5)
hold on
semilogy(freq,abs(FRFm_nHa),'LineWidth',1.5)
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Magnitude of the vertical displacement and acceleration of node H')
legend('Displacement of node H using FE method', 'Acceleration of node H using FE method', 'Displacement of node H using modal approach (first 2 natural freq.)', 'Acceleration of node H using modal approach (first 2 natural freq.)')

subplot(2,1,2)
grid on
plot(freq,angle(FRF_nHh),'LineWidth',1.5)
hold on
grid on
plot(freq,angle(FRF_nHa),'LineWidth',1.5)
hold on
plot(freq,angle(FRFm_nHh),'LineWidth',1.5)
hold on
plot(freq,angle(FRFm_nHa),'LineWidth',1.5)
xlabel('Frequency [Hz]');
ylabel('Phase[rad]');
title('Phase of the vertical displacement and acceleration of node H')
legend('Displacement of node H using FE method', 'Acceleration of node H using FE method', 'Displacement of node H using modal approach (first 2 natural freq.)', 'Acceleration of node H using modal approach (first 2 natural freq.)')



%% 5) Input: vertical displacement A, Output: vertical acceleration F
FY = zeros(ndof,1);
FY(idb(1,2)) = k;
for ii = 1:length(Om)
    X(:,ii) = (-MFF*Om(ii)^2+1i*CFF*Om(ii)+KFF)\FY;
    acc_5(:,ii) = -Om(ii)^2*X(:,ii);
end
FRF_5 = acc_5(idb(6,2),:);

figure(16)
subplot(2,1,1)
semilogy(freq, abs(FRF_5), 'LineWidth',1.5)
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude')
title('Verticale acceleration of node F due to imposed displacement in A')
subplot(2,1,2)
plot(freq, angle(FRF_5),'LineWidth',1.5)
grid on
xlabel('Frequency [Hz]')
ylabel('Phase[rad]')


%% 6) Time story of the steady-state

vel = 12;
lambda1 = 1;
freq_61 = vel/lambda1;
Om_61 = freq_61*2*pi;
amp_61 = 0.001;
F001 = zeros(ndof,1);
F001(idb(1,2)) = k*amp_61;
delay1 = (-0.7-0.42)*2*pi/lambda1; 
F001(idb(2,2)) = k*amp_61*exp(1i*delay1);
X = (-MFF*Om_61^2+1i*CFF*Om_61+KFF)\(F001);
X_H1 = X(idb(8,2));
dt = 0.0001;
t = 0:dt:1;
%a1_t = zeros(length(t),1);
% for ii = 1:length(t)
%     a1_t(ii) = - Om_61^2*abs(X_H)*cos(Om_61*t(ii)+angle(X_H));
% end

% figure(17)
% plot(t,a1_t)
% grid on

lambda2 = 0.6;
freq_62 =vel/lambda2;
Om_62 = freq_62*2*pi;
amp_62 = 0.0005;
F002 = zeros(ndof,1);
F002(idb(1,2)) = k*amp_62;
delay2 = -2*pi*1.12/0.6; 
F002(idb(2,2)) = k*amp_62*exp(1i*delay2);
X = (-MFF*Om_62^2+1i*CFF*Om_62+KFF)\F002;
X_H2 = X(idb(8,2));
dt = 0.0001;
t = 0:dt:1;
%a2_t = zeros(length(t),1);
% for ii = 1:length(t)
% a2_t(ii) =  - Om_62^2*abs(X_H)*cos(Om_62*t(ii)+angle(X_H));
% end

for ii = 1:length(t)
    a_tot(ii) = - Om_61^2 * abs(X_H1).*cos(Om_61*t(ii) + angle(X_H1)) - Om_62^2 * abs(X_H2).*cos(Om_62 *t(ii) + angle(X_H2));  
end

figure(19)
plot(t,a_tot)

%% 7) Static response of the structure
F0 = zeros(ndof,1);
idof_Fsaddle = idb(8,2);
F0(idof_Fsaddle) = -600;
idof_Fhandle = idb(6,2);
F0(idof_Fhandle) = -100;

xF = KFF\F0;

figure(20)
diseg2(xF,100,incid,l,gamma,posit,idb,xy)
title('Static deflection with a scale factor of 100')
