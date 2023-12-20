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
[omega,i_omega] = sort(omega);              % [frord, ordmode] = sort(frq)
freq0 = omega/2/pi;
modes = modes(:,i_omega);

nmodes = 4;
scale_factor = 2;
for ii = 1:nmodes
    mode = modes(:,ii);
    figure(ii+1);
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy)
    legend('Undeformed structure',['Mode ',num2str(ii),' - freq ',num2str(freq0(ii)),' Hz']);
end

%% Damping Matrix

alpha = 6;      % [s^-1]
beta = 1e-5;    % [s]

C = alpha*Mtot + beta*Ktot;
CFF = C(1:ndof, 1:ndof);

%% Frequency Response Function Fc

% Fc related to vertical displacement and acceleration of node F

freq = [0:0.1:200]';
Om = 2*pi*freq;go
F0 = zeros(ndof,1);

idof_j = idb(6,2); % vertical displacement of node F
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end

FRF = X(idof_j,:);
figure(6)
grid on
subplot(2,1,1)
semilogy(freq, abs(FRF),'linewidth',1.5)
grid on
xlabel('frequency [Hz]')
ylabel('abs(FRF)')
subplot(2,1,2)
plot(freq, angle(FRF),'linewidth',1.5)
grid on

% Acceleration of node F
for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
    acc(:,ii) = -Om(ii)^2*X(:,ii);
end

FRF = acc(idof_j,:);
figure(7)
hold on 
grid on
subplot(2,1,1);
semilogy(freq, abs(FRF));
subplot(2,1,2);
plot(freq, angle(FRF));

% Horizontal displacement and acceleration of node H

F0 = zeros(ndof,1);

idof_j = idb(8,1); % horizontal displacement of node H
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end

FRF = X(idof_j,:);
figure(8)
hold on 
grid on
subplot(2,1,1);
semilogy(freq, abs(FRF));
subplot(2,1,2);
plot(freq, angle(FRF));

% Acceleration of node H (horizontal displacement of node H is known)
acc(:,ii) = -Om(ii)^2*X(:,ii);

FRF = acc(idof_j,:);
figure(9)
hold on 
grid on
subplot(2,1,1)
semilogy(freq, abs(FRF));
subplot(2,1,2)
plot(freq, angle(FRF));



