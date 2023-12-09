close all
clear all
clc

%% Parameters

L= 1.2;         %[m]
h= 0.008;       %[m]
b= 0.04;        %[m]
rho= 2700;      %[kg/m^3]
E= 68000000000; %[pa]
J= b*h^3/12;    %[m^4]
m= rho*b*h;     %[kg]

%% Vibration modes

%  standing wave solution
% w=[A*cos(g*x)+B*sin(g*x)+C*cosh(g*x)+D*sinh(g*x)]*cos(omega*t+phi);
% dw/dx = [g*(-A*sin(g*x)+B*cos(g*x)+C*sinh(g*x)+D*cosh(g*x)]*cos(omega*t+phi);
% d^2w/dx^2 = [g^2*(-A*cos(g*x)-B*sin(g*x)+C*cosh(g*x)+D*sinh(g*x)]*cos(omega*t+phi);
% d^3w/dx^3 = [g^3*(A*sin(g*x)-B*cos(g*x)+C*sinh(g*x)+D*cosh(g*x)]*cos(omega*t+phi);


% Boundary conditions 
fs=0.01:0.01:200;
w=2*pi*fs;
g= ((m*w.^2)/(E*J)).^(1/4);

%S= zeros(size(w,2),1);

f=  (-2-2.*cos(g.*L).*cosh(g.*L)).*g.^6 ;




%for i= 1:size(w,2)
 %   H= [1 0 1 0; 0 g(i) 0 g(i); -g(i)^2*cos(g(i)*L) -g(i)^2*sin(g(i)*L) g(i)^2*cosh(g(i)*L) g(i)^2*sinh(g(i)*L); g(i)^3*sin(g(i)*L) -g(i)^3*cos(g(i)*L) g(i)^3*sinh(g(i)*L) g(i)^3*cosh(g(i)*L)];
  %  S(i)= det(H);
%end

%figure;
%semilogy (w,abs(S))

figure;
semilogy (fs,abs(f))
title('Absolute value of det(H) function')

%figure
%plot(w,f)

fs_1= 4.5;
fs_2= 28.2;
fs_3= 79;
fs_4= 154.9;

wn_1= 2*pi*fs_1;
wn_2= 2*pi*fs_2;
wn_3= 2*pi*fs_3;
wn_4= 2*pi*fs_4;

g_1=((m*wn_1.^2)/(E*J)).^(1/4);
g_2=((m*wn_2.^2)/(E*J)).^(1/4);
g_3=((m*wn_3.^2)/(E*J)).^(1/4);
g_4=((m*wn_4.^2)/(E*J)).^(1/4);


wn_vector=[wn_1; wn_2; wn_3; wn_4];
g_vector=[g_1; g_2; g_3; g_4];
constants_matrix= [zeros(4,4)];

for i=1:4 
    H= [1 0 1 0; 0 g_vector(i) 0 g_vector(i); -g_vector(i)^2*cos(g_vector(i)*L) -g_vector(i)^2*sin(g_vector(i)*L) g_vector(i)^2*cosh(g_vector(i)*L) g_vector(i)^2*sinh(g_vector(i)*L); g_vector(i)^3*sin(g_vector(i)*L) -g_vector(i)^3*cos(g_vector(i)*L) g_vector(i)^3*sinh(g_vector(i)*L) g_vector(i)^3*cosh(g_vector(i)*L)];
    N=[H(2,1) H(3,1) H(4,1)]';

    H_hat=[H(2,2) H(2,3) H(2,4); H(3,2) H(3,3) H(3,4); H(4,2) H(4,3) H(4,4)];
    X_hat= -inv(H_hat)*N;
    X= [1; X_hat];
    constants_matrix(:,i)=X;
end

X1= constants_matrix(:,1);
X2= constants_matrix(:,2);
X3= constants_matrix(:,3);
X4= constants_matrix(:,4);

disp(['X1 =' mat2str(X1)]);
disp(['X2 =' mat2str(X2)]);
disp(['X3 =' mat2str(X3)]);
disp(['X4 =' mat2str(X4)]);

x=0:0.01:L;

eigenfunctions= zeros(1,4);
figure;
for i= 1:4;
    eigenfunctions=[constants_matrix(1,i)*cos(g_vector(i)*x)+constants_matrix(2,i)*sin(g_vector(i)*x)+constants_matrix(3,i)*cosh(g_vector(i)*x)+constants_matrix(4,i)*sinh(g_vector(i)*x)];

    subplot(4,1,i);
    plot(x,eigenfunctions/max(abs(eigenfunctions)))
    title('Mode shape')
    xlabel('w')
    ylabel('x')

end


%% Frequency response functions 

%Computation of mi
csi= 0.01;
mass= zeros(1,4);

for i= 1:4 
    eigenfunctions=[constants_matrix(1,i)*cos(g_vector(i)*x)+constants_matrix(2,i)*sin(g_vector(i)*x)+constants_matrix(3,i)*cosh(g_vector(i)*x)+constants_matrix(4,i)*sinh(g_vector(i)*x)];
    Y= m.*eigenfunctions.^2;
    mass(i)= trapz(x,Y);
end

xk=[1.2;1.2;1.2;1.2];
xj=[0.2;0.605;0.8;1.2];                
%xk=[1.2;1;0.605;0.2];
%FRF n°1

eigenfunctions_ktot1= zeros(1,4);
eigenfunctions_jtot1= zeros(1,4);

for l= 1:4 
   
    eigenfunctions_k=-[constants_matrix(1,l)*cos(g_vector(l)*xk(1))+constants_matrix(2,l)*sin(g_vector(l)*xk(1))+constants_matrix(3,l)*cosh(g_vector(l)*xk(1))+constants_matrix(4,l)*sinh(g_vector(l)*xk(1))];
    eigenfunctions_j=[constants_matrix(1,l)*cos(g_vector(l)*xj(1))+constants_matrix(2,l)*sin(g_vector(l)*xj(1))+constants_matrix(3,l)*cosh(g_vector(l)*xj(1))+constants_matrix(4,l)*sinh(g_vector(l)*xj(1))];
    
    eigenfunctions_ktot1(l)= eigenfunctions_k;
    eigenfunctions_jtot1(l)= eigenfunctions_j;
end



Gjk_1= [((eigenfunctions_ktot1(1).*eigenfunctions_jtot1(1))./mass(1))./(-w.^2+(2.*csi.*wn_vector(1).*w).*sqrt(-1)+ wn_vector(1).^2)]+[((eigenfunctions_ktot1(2).*eigenfunctions_jtot1(2))./mass(2))./(-w.^2+(2.*csi.*wn_vector(2).*w).*sqrt(-1)+ wn_vector(2).^2)]+[((eigenfunctions_ktot1(3).*eigenfunctions_jtot1(3))./mass(3))./(-w.^2+(2.*csi.*wn_vector(3).*w).*sqrt(-1)+ wn_vector(3).^2)]+[((eigenfunctions_ktot1(4).*eigenfunctions_jtot1(4))./mass(4))./(-w.^2+(2.*csi.*wn_vector(4).*w).*sqrt(-1)+ wn_vector(4).^2)];


figure;
title('Amplitudes and phases of the FRFs')
subplot(4,2,1)
semilogy(fs,abs(Gjk_1));
title('FRF 1 amplitude')

subplot(4,2,2)
plot(fs,angle(Gjk_1));
title('FRF 1 phase')

%FRF n°2

eigenfunctions_ktot2= zeros(1,4);
eigenfunctions_jtot2= zeros(1,4);

for l= 1:4 
   
    eigenfunctions_k=-[constants_matrix(1,l)*cos(g_vector(l)*xk(2))+constants_matrix(2,l)*sin(g_vector(l)*xk(2))+constants_matrix(3,l)*cosh(g_vector(l)*xk(2))+constants_matrix(4,l)*sinh(g_vector(l)*xk(2))];
    eigenfunctions_j=[constants_matrix(1,l)*cos(g_vector(l)*xj(2))+constants_matrix(2,l)*sin(g_vector(l)*xj(2))+constants_matrix(3,l)*cosh(g_vector(l)*xj(2))+constants_matrix(4,l)*sinh(g_vector(l)*xj(2))];
    
    eigenfunctions_ktot2(l)= eigenfunctions_k;
    eigenfunctions_jtot2(l)= eigenfunctions_j;
end



Gjk_2= [((eigenfunctions_ktot2(1).*eigenfunctions_jtot2(1))./mass(1))./(-w.^2+(2.*csi.*wn_vector(1).*w).*sqrt(-1)+ wn_vector(1).^2)]+[((eigenfunctions_ktot2(2).*eigenfunctions_jtot2(2))./mass(2))./(-w.^2+(2.*csi.*wn_vector(2).*w).*sqrt(-1)+ wn_vector(2).^2)]+[((eigenfunctions_ktot2(3).*eigenfunctions_jtot2(3))./mass(3))./(-w.^2+(2.*csi.*wn_vector(3).*w).*sqrt(-1)+ wn_vector(3).^2)]+[((eigenfunctions_ktot2(4).*eigenfunctions_jtot2(4))./mass(4))./(-w.^2+(2.*csi.*wn_vector(4).*w).*sqrt(-1)+ wn_vector(4).^2)];


subplot(4,2,3)
semilogy(fs,abs(Gjk_2));
title('FRF 2 amplitude')
subplot(4,2,4)
plot(fs,angle(Gjk_2));
title('FRF 2 phase')

%FRF n°3

eigenfunctions_ktot3= zeros(1,4);
eigenfunctions_jtot3= zeros(1,4);

for l= 1:4 
   
    eigenfunctions_k=-[constants_matrix(1,l)*cos(g_vector(l)*xk(3))+constants_matrix(2,l)*sin(g_vector(l)*xk(3))+constants_matrix(3,l)*cosh(g_vector(l)*xk(3))+constants_matrix(4,l)*sinh(g_vector(l)*xk(3))];
    eigenfunctions_j=[constants_matrix(1,l)*cos(g_vector(l)*xj(3))+constants_matrix(2,l)*sin(g_vector(l)*xj(3))+constants_matrix(3,l)*cosh(g_vector(l)*xj(3))+constants_matrix(4,l)*sinh(g_vector(l)*xj(3))];
    
    eigenfunctions_ktot3(l)= eigenfunctions_k;
    eigenfunctions_jtot3(l)= eigenfunctions_j;
end



Gjk_3= [((eigenfunctions_ktot3(1).*eigenfunctions_jtot3(1))./mass(1))./(-w.^2+(2.*csi.*wn_vector(1).*w).*sqrt(-1)+ wn_vector(1).^2)]+[((eigenfunctions_ktot3(2).*eigenfunctions_jtot3(2))./mass(2))./(-w.^2+(2.*csi.*wn_vector(2).*w).*sqrt(-1)+ wn_vector(2).^2)]+[((eigenfunctions_ktot3(3).*eigenfunctions_jtot3(3))./mass(3))./(-w.^2+(2.*csi.*wn_vector(3).*w).*sqrt(-1)+ wn_vector(3).^2)]+[((eigenfunctions_ktot3(4).*eigenfunctions_jtot3(4))./mass(4))./(-w.^2+(2.*csi.*wn_vector(4).*w).*sqrt(-1)+ wn_vector(4).^2)];


subplot(4,2,5)
semilogy(fs,abs(Gjk_3));
title('FRF 3 amplitude')
subplot(4,2,6)
plot(fs,angle(Gjk_3));
title('FRF 3 phase')


%FRF n°4

eigenfunctions_ktot4= zeros(1,4);
eigenfunctions_jtot4= zeros(1,4);

for l= 1:4 
   
    eigenfunctions_k=-[constants_matrix(1,l)*cos(g_vector(l)*xk(4))+constants_matrix(2,l)*sin(g_vector(l)*xk(4))+constants_matrix(3,l)*cosh(g_vector(l)*xk(4))+constants_matrix(4,l)*sinh(g_vector(l)*xk(4))];
    eigenfunctions_j=[constants_matrix(1,l)*cos(g_vector(l)*xj(4))+constants_matrix(2,l)*sin(g_vector(l)*xj(4))+constants_matrix(3,l)*cosh(g_vector(l)*xj(4))+constants_matrix(4,l)*sinh(g_vector(l)*xj(4))];
    
    eigenfunctions_ktot4(l)= eigenfunctions_k;
    eigenfunctions_jtot4(l)= eigenfunctions_j;
end



Gjk_4= [((eigenfunctions_ktot4(1).*eigenfunctions_jtot4(1))./mass(1))./(-w.^2+(2.*csi.*wn_vector(1).*w).*sqrt(-1)+ wn_vector(1).^2)]+[((eigenfunctions_ktot4(2).*eigenfunctions_jtot4(2))./mass(2))./(-w.^2+(2.*csi.*wn_vector(2).*w).*sqrt(-1)+ wn_vector(2).^2)]+[((eigenfunctions_ktot4(3).*eigenfunctions_jtot4(3))./mass(3))./(-w.^2+(2.*csi.*wn_vector(3).*w).*sqrt(-1)+ wn_vector(3).^2)]+[((eigenfunctions_ktot4(4).*eigenfunctions_jtot4(4))./mass(4))./(-w.^2+(2.*csi.*wn_vector(4).*w).*sqrt(-1)+ wn_vector(4).^2)];


subplot(4,2,7)
semilogy(fs,abs(Gjk_4));
title('FRF 4 amplitude')
subplot(4,2,8)
plot(fs,angle(Gjk_4));
title('FRF 4 phase')

%Comparison between FRFs

figure
semilogy(fs,abs(Gjk_1));
hold on
semilogy(fs,abs(Gjk_2));
hold on
semilogy(fs,abs(Gjk_3));
hold on
semilogy(fs,abs(Gjk_4));

title('comparison between FRFs amplitudes')
legend('FRF1','FRF2','FRF3','FRF4')

% figure
% plot(fs,angle(Gjk_1));
% hold on
% plot(fs,angle(Gjk_2));
% hold on
% plot(fs,angle(Gjk_3));
% hold on
% plot(fs,angle(Gjk_4));
% 
% title('comparison between FRFs phases')
% legend('FRF1','FRF2','FRF3','FRF4')
Gjk_tot=[Gjk_1;Gjk_2;Gjk_3;Gjk_4];
fs_tot=[fs_1;fs_2;fs_3;fs_4];
w2= 4.545*2*pi;
w1= 4.455*2*pi;
est_damp1= (w2^2-w1^2)/(4*wn_1^2);

w2= 79.74*2*pi;
w1= 78.17*2*pi;
est_damp3= (w2^2-w1^2)/(4*wn_3^2);

w2= 28.46*2*pi;
w1= 27.88*2*pi;
est_damp2= (w2^2-w1^2)/(4*wn_2^2);

w2= 156.36*2*pi;
w1= 153.31*2*pi;
est_damp4= (w2^2-w1^2)/(4*wn_4^2);

est_damp= [est_damp1;est_damp2;est_damp3;est_damp4];

entrata_Gjk= fs_tot(1)/0.01;
estim11= -imag(Gjk_1(entrata_Gjk)).*2.*wn_1.^2.*est_damp1;
estim12= -imag(Gjk_2(entrata_Gjk)).*2.*wn_1.^2.*est_damp1;
estim13= -imag(Gjk_3(entrata_Gjk)).*2.*wn_1.^2.*est_damp1;
estim14= -imag(Gjk_4(entrata_Gjk)).*2.*wn_1.^2.*est_damp1;

entrata_Gjk= fs_tot(2)/0.01;
estim21= -imag(Gjk_1(entrata_Gjk)).*2.*wn_2.^2.*est_damp2;
estim22= -imag(Gjk_2(entrata_Gjk)).*2.*wn_2.^2.*est_damp2;
estim23= -imag(Gjk_3(entrata_Gjk)).*2.*wn_2.^2.*est_damp2;
estim24= -imag(Gjk_4(entrata_Gjk)).*2.*wn_2.^2.*est_damp2;

entrata_Gjk= fs_tot(3)/0.01;
estim31= -imag(Gjk_1(entrata_Gjk)).*2.*wn_3.^2.*est_damp3;
estim32= -imag(Gjk_2(entrata_Gjk)).*2.*wn_3.^2.*est_damp3;
estim33= -imag(Gjk_3(entrata_Gjk)).*2.*wn_3.^2.*0;
estim34= -imag(Gjk_4(entrata_Gjk)).*2.*wn_3.^2.*0;

entrata_Gjk= fs_tot(4)/0.01;
estim41= -imag(Gjk_1(entrata_Gjk)).*2.*wn_4.^2.*est_damp4;
estim42= -imag(Gjk_2(entrata_Gjk)).*2.*wn_4.^2.*est_damp4;
estim43= -imag(Gjk_3(entrata_Gjk)).*2.*wn_4.^2.*est_damp4;
estim44= -imag(Gjk_4(entrata_Gjk)).*2.*wn_4.^2.*est_damp4;

estimation_A= [estim11 estim21 estim31 estim41;     %first number=peak, second number=FRF
               estim12 estim22 estim32 estim42;
               estim13 estim23 estim33 estim43;
               estim14 estim24 estim34 estim44];


%% Modal parameters identification
% Peak 1 FRF 1
mode_vector=zeros(4,4);
for k=1:4   %number of peaks

fs=fs_tot(k)-2:0.01:fs_tot(k)+2;
omegafs = 2*pi*fs;

GrEXP = [ Gjk_1((fs_tot(k)-2)*100:(fs_tot(k)-2)*100-1+size(fs,2)) ; Gjk_2((fs_tot(k)-2)*100:(fs_tot(k)-2)*100-1+size(fs,2)) ; Gjk_3((fs_tot(k)-2)*100:(fs_tot(k)-2)*100-1+size(fs,2)) ; Gjk_4((fs_tot(k)-2)*100:(fs_tot(k)-2)*100-1+size(fs,2)) ];

% Least square minimization (have a look on epsilon function)

% Initial guesses
x0 = [estimation_A(1,k),estimation_A(2,k),estimation_A(3,k),estimation_A(4,k),est_damp(k),wn_vector(k),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% Lsqm function
%options=optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',13000);
x = lsqnonlin(@(x) epsilon1(x, omegafs, GrEXP), x0,[],[],[]);
% Having a first look at x
disp(x);

% Computing of GrNUMs to compare them to Gjk

%  each GrNUMi is supposed to match with the peak number i 
for s=1:4   %number of FRFs
GrNUM1 = (x(s)./(-omegafs.^2 + 1i*2*x(5)*x(6)*omegafs + x(6)^2)) + (x(s+6)+1i*x(s+10))./(omegafs.^2) + x(s+14)+1i*x(s+18);

mode_entr= imag(GrNUM1(1,201));
mode_vector(k,s)= mode_entr;                      %raw=peak, colomn=frf

freq_range= 0.01:0.01:200;
figure
subplot(2,1,1)
semilogy(freq_range,abs(Gjk_tot(s,:)))
hold on 
semilogy(fs,abs(GrNUM1),'o')
subplot(2,1,2)
plot(freq_range,angle(Gjk_tot(s,:)))
hold on 
plot(fs,angle(GrNUM1),'o')
end
end


% 
%  fs=115:0.01:200;
% 
% y_damp= max(abs(Gjk_1(11500:11499+size(fs,2))))*sqrt(2)/2;
% figure
% semilogy(fs,abs(Gjk_1(11500:11499+size(fs,2))));
% hold on
% k=y_damp*ones(size(fs,2));
% plot(fs,k);

%Mode shape identification
X1num= mode_vector(1,:)*mass(1);
figure

x=0:0.01:L;
eigenfunctions=[constants_matrix(1,1)*cos(g_vector(1)*x)+constants_matrix(2,1)*sin(g_vector(1)*x)+constants_matrix(3,1)*cosh(g_vector(1)*x)+constants_matrix(4,1)*sinh(g_vector(1)*x)];
plot(xj,X1num,'o')
hold on
plot(x,eigenfunctions)
hold on


% figure 
% plot(xj,x2(1:4),'o')
% 
% figure 
% plot(xj,x3(1:4),'o')
% 
% figure 

% 
% figure 
% p
