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
% dis_stru(posit,l,gamma,xy,pr,idb,ndof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble mass and stiffness matricies
[M,K] = assem(incid,l,m,EA,EJ,T,gamma,idb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add concentrated elements
% Modal matrices

% Concentrated spring in node 1
k = 4e+05;
i_ndof_spring = idb(1,2);
E_k1 = zeros(1,45);
E_k1(:, i_ndof_spring) = 1;
Kx1 = E_k1'*k*E_k1;

% Concentrated spring in node 2
k = 4e+05;
i_ndof_spring = idb(2,2);
E_k1 = zeros(1,45);
E_k1(:, i_ndof_spring) = 1;
Kx2 = E_k1'*k*E_k1;

% Ktot
Ktot = K + Kx1 + Kx2;

% Concentrated masses
m1 = 1; %[kg]
m2 = 2; %[kg]
m3 = 2; %[kg]
m4 = 1.5; %[kg]

% m1
M1_h = diag([m1 m1 0]);
i_ndof_m1 = idb(1,:);
E_m1 = zeros(3,45);
E_m1(:,i_ndof_m1) = eye(3);
M1 = E_m1'*M1_h*E_m1;

% m2
M2_h = diag([m2 m2 0]);
i_ndof_m2 = idb(3,:);
E_m2 = zeros(3,45);
E_m2(:,i_ndof_m2) = eye(3);
M2 = E_m2'*M2_h*E_m2;

% m3
M3_h = diag([m3 m3 0]);
i_ndof_m3 = idb(2,:);
E_m3 = zeros(3,45);
E_m3(:,i_ndof_m3) = eye(3);
M3 = E_m3'*M3_h*E_m3;

% m4
M4_h = diag([m4 m4 0]);
i_ndof_m4 = idb(6,:);
E_m4 = zeros(3,45);
E_m4(:,i_ndof_m4) = eye(3);
M4 = E_m4'*M4_h*E_m4;

% Mtot
Mtot = M + M1 + M2 + M3 + M4;

%% Compute natural frequencies and mode shapes
MFF = Mtot(1:ndof,1:ndof);
MCF = Mtot(ndof+1:end,1:ndof);
MFC = Mtot(1:ndof,ndof+1:end);
MCC = Mtot(ndof+1:end,ndof+1:end);

KFF = Ktot(1:ndof,1:ndof);
KCF = Ktot(ndof+1:end,1:ndof);
KFC = Ktot(1:ndof,ndof+1:end);
KCC = Ktot(ndof+1:end,ndof+1:end);

[modes, omega] = eig(inv(MFF)*KFF);
omega = sqrt(diag(omega));
[omega,i_omega] = sort(omega);
freq0 = omega/2/pi;
modes = modes(:,i_omega);

% nmodes = 4;
% scale_factor = 2;
% for ii = 1:nmodes
%     mode = modes(:,ii);
%     figure(ii+1);
%     diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy)
%     legend('Undeformed structure',['Mode ',num2str(ii),' - freq ',num2str(freq0(ii)),' Hz']);
% end

%% Damping Matrix

alpha = 6;      % [s^-1]
beta = 1e-5;    % [s]

C = alpha*Mtot + beta*Ktot;
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
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end
FRF_nFv = X(idof_j,:);

% figure(6)
% subplot(2,1,1)
% semilogy(freq, abs(FRF_nFv),'linewidth',1.5)
% grid on
% xlabel('frequency [Hz]')
% ylabel('abs(FRF)')
% title('Vertical displacement of node F')
% subplot(2,1,2)
% plot(freq, angle(FRF_nFv),'linewidth',1.5)
% grid on
% xlabel('frequency [Hz]')
% ylabel('\phi(FRF)')

% Acceleration of node F
for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
    acc(:,ii) = -Om(ii)^2*X(:,ii);
end
FRF_nFa = acc(idof_j,:);

% figure(7)
% subplot(2,1,1)
% grid on
% semilogy(freq, abs(FRF_nFa));
% subplot(2,1,2)
% plot(freq, angle(FRF_nFa));

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

% figure(8)
% hold on 
% grid on
% subplot(2,1,1);
% semilogy(freq, abs(FRF_nHh));
% subplot(2,1,2);
% plot(freq, angle(FRF_nHh));

% Acceleration of node H (horizontal displacement of node H is known)
acc(:,ii) = -Om(ii)^2*X(:,ii);
FRF_nHa = acc(idof_j,:);

% figure(9)
% subplot(2,1,1)
% semilogy(freq, abs(FRF_nHa));
% grid on
% subplot(2,1,2)
% plot(freq, angle(FRF_nHa));

%% 3B) Internal forces in node 11
nel = 7;
F0 = zeros(ndof,1);

idof_7 = idb(7,:); 
idof_11 = idb(11,:);
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end
X_7 = X(idof_7,:);
X_11 = X(idof_11,:);

lambda = [cos(gamma(nel)) sin(gamma(nel)) 0;
        -sin(gamma(nel)) cos(gamma(nel)) 0;
         0                0               1];

L_el = l(nel);
csi = L_el;
loc_X_7 = lambda*X_7; %Xi
loc_X_11 = lambda*X_11; %Xj

% b for axial force N
b = (loc_X_11(1,:) - loc_X_7(1,:))/L_el;

N_11 = EA(nel)*b;

% c and d for shear force T and bending moment M
c = -3/L_el^2*loc_X_7(2,:) + 3/L_el^2*loc_X_11(2,:) - 2/L_el*loc_X_7(3,:) - 1/L_el*loc_X_11(3,:);
d = 2/L_el^3*loc_X_7(2,:) - 2/L_el^3*loc_X_11(2,:) + 1/L_el^2*loc_X_7(3,:) + 1/L_el^2*loc_X_11(3,:);

T_11 = EJ(nel)*(6*d);
M_11 = EJ(nel)*(2*c + 6*d*csi);

% figure(10)
% subplot(2,1,1)
% semilogy(freq,abs(N_11))
% hold on
% semilogy(freq,abs(T_11))
% semilogy(freq,abs(M_11))
% grid on
% subplot(2,1,2)
% plot(freq,angle(N_11))
% hold on
% grid on
% plot(freq,angle(T_11))
% plot(freq,angle(M_11))

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

% figure(11)
% subplot(2,1,1)
% semilogy(freq, abs(FRF_RV));
% grid on
% subplot(2,1,2)
% plot(freq, angle(FRF_RV));
% grid on

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

% figure(12)
% subplot(2,1,1)
% semilogy(freq,abs(FRF_nFv))
% grid on
% hold on
% semilogy(freq,abs(FRFm_nFv))
% subplot(2,1,2)
% plot(freq,angle(FRF_nFv))
% grid on
% hold on
% plot(freq,angle(FRFm_nFv))
% 
% figure(13)
% subplot(2,1,1)
% semilogy(freq,abs(FRF_nFa))
% grid on
% hold on
% semilogy(freq,abs(FRFm_nFa))
% subplot(2,1,2)
% plot(freq,angle(FRF_nFa))
% grid on
% hold on
% plot(freq,angle(FRFm_nFa))

%% 5) Input: vertical displacement A, Output: vertical acceleration F
FY = zeros(ndof,1);
FY(idb(1,2)) = k;
for ii = 1:length(Om)
    % F_temp = (-MFF*Om(ii)^2+1i*CFF*Om(ii)+KFF)*Y;
    X(:,ii) = (-MFF*Om(ii)^2+1i*CFF*Om(ii)+KFF)\FY;
    acc_5(:,ii) = -Om(ii)^2*X(:,ii);
end
FRF_5 = acc_5(idb(6,2),:);

% figure(14)
% subplot(2,1,1)
% semilogy(freq, abs(FRF_5));
% grid on
% subplot(2,1,2)
% plot(freq, angle(FRF_5))
% grid on

%% 6) Time story of the steady-state

vel = 12;
lambda1 = 1;
freq_61 = vel/lambda1;
Om_61 = freq_61*2*pi;
amp_61 = 0.001;
F00 = zeros(ndof,1);
F00(idb(1,2)) = k*amp_61;
delay1 = (-0.7-0.42)*2*pi/lambda1; 
F00(idb(2,2)) = k*amp_61*exp(1i*delay1);
X = (-MFF*Om_61^2+1i*CFF*Om_61+KFF)\F00;
X_H = X(idb(8,2));
dt = 0.0001;
t = 0:dt:1;
a1_t = zeros(length(t),1);
for ii = 1:length(t)
a1_t(ii) = Om_61^2*abs(X_H)*cos(Om_61*t(ii)+angle(X_H));
end

figure()
plot(t,a1_t)

lambda2 = 0.6;
freq_62 =vel/lambda2;
Om_62 = freq_62*2*pi;
amp_62 = 0.0005;
F00 = zeros(ndof,1);
F00(idb(1,2)) = k*amp_62;
delay2 = -2*pi*1.12/0.6; 
F00(idb(2,2)) = k*amp_62*exp(1i*delay2);
X = (-MFF*Om_62^2+1i*CFF*Om_62+KFF)\F00;
X_H = X(idb(8,2));
dt = 0.0001;
t = 0:dt:1;
a2_t = zeros(length(t),1);
for ii = 1:length(t)
a2_t(ii) = -Om_62^2*abs(X_H)*cos(Om_62*t(ii)+angle(X_H));
end
figure()
plot(t,a2_t)

a_tot = a1_t + a2_t;
figure()
plot(t,a_tot)

%% 7) Static response of the structure
F0 = zeros(ndof,1);
idof_Fsaddle = idb(8,2);
F0(idof_Fsaddle) = -600;
idof_Fhandle = idb(6,2);
F0(idof_Fhandle) = -100;

xF = KFF\F0;

% figure()
% diseg2(xF,20,incid,l,gamma,posit,idb,xy)
