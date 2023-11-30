%% Initialisation of the problem

clear all

[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure;

dis_stru(posit,l,gamma,xy,pr,idb,ndof); % for drawing the bike

omegamax = 200*2*pi;
omega1 = ((pi./l).^2).*sqrt(EJ./m);

Verif = omega1/omegamax;

[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

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

figure(1);
scale_factor = 1.5;
for ii = 1:2
    mode = modes_sorted(:,ii);
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy);
end

%% Damping matrix 

alpha = 6;
beta = 10e-5;
CFF = alpha*MFF + beta*KFF;

%% Frequency Response Function

% Get the value in ibd of C ( 3rd node so 3rd row ) for Y displacement (
% 2nd column)

F0 = zeros(ndof,1);
index = idb(3,2);
F0(index) = 1;
om = (0:1:200)*2*pi;
for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

figure(2);
plot(0:1:200,X);









