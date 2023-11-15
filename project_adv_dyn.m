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
fs=0:0.1:200;
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

figure;
for i= 1:4;
    eigenfunctions=[constants_matrix(1,i)*cos(g_vector(i)*x)+constants_matrix(2,i)*sin(g_vector(i)*x)+constants_matrix(3,i)*cosh(g_vector(i)*x)+constants_matrix(4,i)*sinh(g_vector(i)*x)];

    subplot(4,1,i);
    plot(x,eigenfunctions)
    title('Mode shape')
    xlabel('w')
    ylabel('x')

end


%% Frequency response functions 

csi= 0.01;
mass= zeros(1,4);

for i= 1:4 
    eigenfunctions=[constants_matrix(1,i)*cos(g_vector(i)*x)+constants_matrix(2,i)*sin(g_vector(i)*x)+constants_matrix(3,i)*cosh(g_vector(i)*x)+constants_matrix(4,i)*sinh(g_vector(i)*x)];
    Y= m.*eigenfunctions.^2;
    mass(i)= trapz(x,Y);
end

xk= L;
xj= 0.5*L;

eigenfunctions_ktot= zeros(1,4);
eigenfunctions_jtot= zeros(1,4);

for l= 1:4 
    eigenfunctions_k=[constants_matrix(1,l)*cos(g_vector(l)*xk)+constants_matrix(2,l)*sin(g_vector(l)*xk)+constants_matrix(3,l)*cosh(g_vector(l)*xk)+constants_matrix(4,l)*sinh(g_vector(l)*xk)];
    eigenfunctions_j=[constants_matrix(1,l)*cos(g_vector(l)*xj)+constants_matrix(2,l)*sin(g_vector(l)*xj)+constants_matrix(3,l)*cosh(g_vector(l)*xj)+constants_matrix(4,l)*sinh(g_vector(l)*xj)];
    
    eigenfunctions_ktot(l)= eigenfunctions_k;
    eigenfunctions_jtot(l)= eigenfunctions_j;
end


OMEGA= 0:0.1:1000;

Gjk_1= [((eigenfunctions_ktot(1).*eigenfunctions_jtot(1))./mass(1))./(-OMEGA.^2+(2.*csi.*wn_vector(1).*OMEGA).*sqrt(-1)+ wn_vector(1))]+[((eigenfunctions_ktot(2).*eigenfunctions_jtot(2))./mass(2))./(-OMEGA.^2+(2.*csi.*wn_vector(2).*OMEGA).*sqrt(-1)+ wn_vector(2))]+[((eigenfunctions_ktot(3).*eigenfunctions_jtot(3))./mass(3))./(-OMEGA.^2+(2.*csi.*wn_vector(3).*OMEGA).*sqrt(-1)+ wn_vector(3))]+[((eigenfunctions_ktot(4).*eigenfunctions_jtot(4))./mass(4))./(-OMEGA.^2+(2.*csi.*wn_vector(4).*OMEGA).*sqrt(-1)+ wn_vector(4))];


figure;
plot(OMEGA,Gjk_1);


  










   


