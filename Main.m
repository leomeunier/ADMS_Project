%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical System Dynamics
% FEM script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
clc;
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
[omega,i_omega] = sort(omega);              % [frord, ordmode] = sort(frq)
freq0 = omega/2/pi;
modes = modes(:,i_omega);

nmodes = 4;
scale_factor = 2;
for ii=1:nmodes
    mode = modes(:,ii);
    figure();
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy)
    title(['mode ',num2str(ii),' freq ',num2str(freq0(ii)),' Hz']);
end

% Damping Matrix

alpha = 6;      % [s^-1]
beta = 1e-5;    % [s]

C = alpha*M + beta*K;
CFF = C(1:ndof, 1:ndof);

%% Frequency Response Function Fc

% Fc related to vertical displacement and acceleration of node F

freq = [0:0.1:200]';
Om = 2*pi*freq;
F0 = zeros(ndof,1);

idof_j = idb(6,2); % vertical displacement of node F
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end

FRF = X(idof_j,:);
figure()
hold on 
grid on
subplot(2,2,1);
semilogy(freq, abs(FRF));
subplot(2,1,2);
plot(freq, angle(FRF));

% Manca accelerazione


% Horizontal displacement and acceleration of node H

F0 = zeros(ndof,1);

idof_j = idb(8,2); % vertical displacement of node H
idof_k = idb(3,2); % force in node C

F0(idof_k) = 1;

for ii = 1:length(Om)
    A = -Om(ii)^2*MFF + 1i*Om(ii)*CFF + KFF;
    X(:,ii) = A\F0;      % Matrix in which each row is a FRF of dof
end

FRF = X(idof_j,:);
figure()
hold on 
grid on
subplot(2,2,1);
semilogy(freq, abs(FRF));
subplot(2,1,2);
plot(freq, angle(FRF));

% Manca accelerazione





















