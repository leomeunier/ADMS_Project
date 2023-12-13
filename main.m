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

[modes, omega2] = eig(inv(MFF)*KFF);
omega = diag(sqrt(omega2));
% Sort frequencies in ascending order
[omega_sorted, omega_sorted_indices] = sort(omega);
% Sort mode shapes in ascending order
modes_sorted = modes(:,omega_sorted_indices);

disp(omega_sorted/(2*pi));

%% Drawing modes shapes

figure(2);
scale_factor = 1.5;
for ii = 1:4
    mode = modes_sorted(:,ii);
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy);
end

%% Damping matrix 

alpha = 6;
beta = 10e-5;
CFF = alpha*MFF + beta*KFF;

%% Frequency Response Function

% Get the value in idb of C ( 3rd node so 3rd row ) for Y displacement (
% 2nd column)

F0 = zeros(ndof,1);
index = idb(3,2);
F0(index) = 1;
fs = 0:0.1:200;
om = (0:0.1:200)*2*pi;
for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

% Vertical displacement of node F

vf = idb(10,2);

figure(3)
subplot(2,1,1)
semilogy(fs,abs(X(vf,:)'))
subplot(2,1,2)
plot(fs,angle(X(vf,:)'))

%  Horizontal displacement of node H

hh = idb(13,1);

figure(5)
subplot(2,1,1)
semilogy(fs,abs(X(hh,:)'))
subplot(2,1,2)
plot(fs,angle(X(hh,:)'))







