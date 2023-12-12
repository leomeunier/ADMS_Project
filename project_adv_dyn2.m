close all
clear all
clc
%properties
red = [20e-3 1e-3];
blue = [32e-3 1.5e-3];
green = [50e-3 2.5e-3];

Em = 7e10;
rho = 2700;

colors = [ red ; blue ; green];
m = zeros(1,3);
EA = zeros(1,3);
EJ = zeros(1,3);
i = 0;

for i = 1:3
    m(i)= rho*pi*(colors(i,1)^2 - (colors(i,1)-colors(i,2))^2)/4;
    EA(i)= Em*pi*(colors(i,1)^2 - (colors(i,1)-colors(i,2))^2)/4;
    EJ(i)= Em*(pi/64)*((colors(i,1))^4 - (colors(i,1)-colors(i,2))^4);
end

c2 = min(EJ./m);
c = sqrt(c2);
fmax = 200;
Lmax = sqrt((pi^2/(2.5*2*pi*fmax))*c);

A = [0.7 0];
B = [-0.42 0];
C = [0 0];
D = [0.56 0.42];
E = [0.52 0.54];
F = [0.49 0.63];
G = [-0.08 0.48];
H = [-0.12 0.72];

%%

[ADMS_FEM_01,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure;

[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

dis_stru(posit,l,gamma,xy,pr,idb,ndof);

L_elements=[norm(xy(1,:)-xy(2,:)),norm(xy(2,:)-xy(3,:)),norm(xy(3,:)-xy(4,:)),norm(xy(4,:)-xy(5,:)),norm(xy(5,:)-xy(6,:)),norm(xy(7,:)-xy(6,:)),norm(xy(8,:)-xy(7,:)),norm(xy(6,:)-xy(9,:)),norm(xy(9,:)-xy(10,:)),norm(xy(9,:)-xy(11,:)),norm(xy(11,:)-xy(12,:)),norm(xy(12,:)-xy(13,:)),norm(xy(12,:)-xy(14,:)),norm(xy(12,:)-xy(15,:)),norm(xy(14,:)-xy(1,:)),norm(xy(15,:)-xy(3,:))];

Wn_elements=zeros(1,16);

red_vec=[1,2,13,15];
blue_vec=[6:12,14,16];
green_vec=[3:5];

for ii=1:4
    W=(pi/L_elements(red_vec(ii)))^2*sqrt(EJ(1)/m(1));
    Wn_elements(1,red_vec(ii))=W;
end

for ii=1:9
    W=(pi/L_elements(blue_vec(ii)))^2*sqrt(EJ(6)/m(6));
    Wn_elements(1,blue_vec(ii))=W;
end


for ii=1:3
    W=(pi/L_elements(green_vec(ii)))^2*sqrt(EJ(3)/m(3));
    Wn_elements(1,green_vec(ii))=W;
end

ratios=Wn_elements./(fmax*2*pi);

if ratios>2
    disp(['the elements are in quasi static region'])
else
    disp(['the elements are not in quasi static region'])
end

MFF=M(1:ndof,1:ndof);
MFC=M(1:ndof,ndof+1:45);
MCF=M(ndof+1:45,1:ndof);
MCC=M(ndof+1:45,ndof+1:45);

KFF=K(1:ndof,1:ndof);
KFC=K(1:ndof,ndof+1:45);
KCF=K(ndof+1:45,1:ndof);
KCC=K(ndof+1:45,ndof+1:45);

[modes omega2] = eig(inv(MFF)*KFF);
omega = diag(sqrt(omega2));
% Sort frequencies in ascending order
[omega_sorted omega_sorted_indices] = sort(omega);
% Sort mode shapes in ascending order
modes_sorted = modes(:,omega_sorted_indices);

Wn_structure= omega_sorted(2:5);

scale_factor=1.5;
for ii=2:5
    mode=modes_sorted(:,ii);
    figure
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy);
end

alpha=6;
beta=1e-5;
CFF=alpha*MFF+beta*KFF;

F0 = zeros(ndof,1);
index = idb(3,2);
F0(index) = 1;
om = (0:0.1:200)*2*pi;
for ii=1:length(om)
A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
X(:,ii) = A\F0;
end

displ_F= idb(10,2);
X_f= X(displ_F,:);
figure
subplot(2,1,1)
semilogy((0:0.1:200), abs(X_f))
subplot(2,1,2)
plot((0:0.1:200),angle(X_f))

