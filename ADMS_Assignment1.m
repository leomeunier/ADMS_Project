clc
clear
close all

% ADMS_ASSIGNMENT1_PartA

%% Data

L = 1.2;         % [m] Length
h = 0.008;       % [m] Thickness
b = 0.04;        % [m] Width
rho = 2700;      % [km/m3] Density
E = 68e+9;       % [Pa] Young's Modulus
m = rho*b*h;     % [kg/m] mass
J = (b*h^3)/12;  % [m^4]

%% Matrix formulation

f = linspace(0, 200, 1000000);
omega = 2*pi*f;
gamma = ((m*(omega.^2))/(E*J)).^(1/4);

H = @(x) [          1            0              1                0;
                    0            1              0                1;
               sin(x*L)      -cos(x*L)      sinh(x*L)        cosh(x*L);
               -cos(x*L)     -sin(x*L)      cosh(x*L)        sinh(x*L)] ;

dets =[];

for i=1:length(gamma)
    dets(i) = det(H(gamma(i)));
end

semilogy(f,abs(dets), 'LineWidth',2)
hold on 
grid on

gamma0 = [];

for i = 2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        gamma0 = [gamma0 , gamma(i)];
    end
end

omega0 = [];

for i = 1:length(gamma0)
    omega0(i) = sqrt((gamma0(i)^4)*(E*J)/m); % Natural pulsations
end

freq = [];

for i = 1:4
    freq(i) = omega0(i)/(2*pi); % Natural frequencies: frequencies at which my system vibrates with a certain shape
end


%% Mode shapes = shape that my sysem assumes when it vibrates at a certain natural frequency, it's a picture of the system

Hgamma = [];
z = ones(4,4);

for i = 1:length(omega0)
    Hgamma = H(gamma0(i));
    H_new = Hgamma(2:end, 2:end);
    N = Hgamma(2:end,1);
    k = H_new\(-N);
    z(2:end, i) = k; % I found the coefficients needed to evaluate the transverse vibration w.
    % z is a matrix in which each column is a vector with the coefficients
    % for each natural frequency.
end

for i = 1:4
    phi = @(x) z(1,i)*cos(gamma0(i)*x) + z(2,i)*sin(gamma0(i)*x) + z(3,i)*cosh(gamma0(i)*x) + z(4,i)*sinh(gamma0(i)*x);
    tratt = @(x) x*0;
    distance = linspace(0,L,1000);
    maxi(i) = (max(abs(phi(distance))));
    phinorm = phi(distance)/maxi(i);
    PHInorm(i,:) = phinorm;

    if i == 1
        phi1 = phi;
    end

    if i == 2
        phi2 = phi;
    end

    if i == 3
        phi3 = phi;
    end

    if i == 4
        phi4 = phi;
    end

    % I obtained the mode shapes corresponding to the first 4 natural
    % frequencies of the system.

    figure(i+1)
    plot(distance, phinorm, 'LineWidth',2)
    hold on, grid on
    axis([0 L  -1 1])
    plot(distance,tratt(distance),'k--')
    xlabel('Distance x from the fixed end [m]')
    ylabel('Mode shapes \phi (x)'); 
    freqnow = freq(i);
    legend(sprintf('mode %g f_{0} = %g', i, freqnow));
end


%% Frequency response function

% Positions to compare my FRF to the one on slide 7

x_j1 = 0.2;            %[m] output position
x_k1 = 1.2;            %[m] input position
csi = ones(4,1)/100;   % damping ratio

% Modal mass
xx = linspace(0, L, 1000);

for ii = 1:length(omega0)
    mm(ii) = trapz(xx, m*PHInorm(ii, :).^2);
end

G_jk = @(OM,x_j,x_k)   (-phi1(x_j)'.*phi1(x_k)./((maxi(1)^2*mm(1)).*(-OM.^2 + 1i*2*csi(1)*omega0(1)*OM+omega0(1)^2))+...
                        -phi2(x_j)'.*phi2(x_k)./((maxi(2)^2*mm(2)).*(-OM.^2 + 1i*2*csi(2)*omega0(2)*OM+omega0(2)^2))+...
                        -phi3(x_j)'.*phi3(x_k)./((maxi(3)^2*mm(3)).*(-OM.^2 + 1i*2*csi(3)*omega0(3)*OM+omega0(3)^2))+...
                        -phi4(x_j)'.*phi4(x_k)./((maxi(4)^2*mm(4)).*(-OM.^2 + 1i*2*csi(4)*omega0(4)*OM+omega0(4)^2)));

% I transpose phi so that I can take many outputs at a time, the result
% will be r rows corresponding to r FRF and r outputs.

OM_vect = 2*pi*(0:0.02:200);

figure (6)
subplot(2,1,1)
semilogy(OM_vect/2/pi, abs(G_jk(OM_vect,x_j1,x_k1)), 'LineWidth', 2)
ylim([1e-6 0.3e-1])
xlabel('Frequency [Hz]')
ylabel('[m/N]'); 
legend('G_{jk} ( \Omega ) = x_{j}/F_{k}, x_{j} = 0.2 m , x_{k} = 1.2 m')
grid on; hold on;

subplot(2,1,2)
plot(OM_vect/2/pi, angle(G_jk(OM_vect,x_j1,x_k1)), 'LineWidth',2)
grid on

% TO DO: 
% Other combinations of input-output positions
% Suggestions by La Paglia: 
% colocated FRF (input = output)
% observability: place the sensor in a vibration node (x_j)
% controllability: place the force in a vibration node (x_k)
% reciprocity

% Altre G_jk
xk=[1.2;0.2;0.605;1];
xj=[0.2;1.2;0.8;0.605];


for i = 2:length(xj)
    figure(5+i)
    subplot(2,1,1)
    semilogy(OM_vect/2/pi, abs(G_jk(OM_vect,xj(i),xk(i))), 'LineWidth', 2)
    xlabel('Frequency [Hz]')
    ylabel('[m/N]'); 
   % legend('G_{jk} ( \Omega ) = x_{j}/F_{k}, x_{j} = %g , x_{k} = %g', xj(i),xk(i))
    grid on; hold on;
    subplot(2,1,2)
    plot(OM_vect/2/pi, angle(G_jk(OM_vect,xj(i),xk(i))), 'LineWidth',2)
    grid on
    hold on
end



%% Numerical Frequency Response Function - Comparison with the one on slide 12

% MODE 1 

x_k=[1.2;0.2;0.605;1];
x_j=[0.2;1.2;0.8;0.605];

for s = 1:length(x_j)

G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;

OM_vect = 2*pi*(4:0.005:5); 

err = @ (X) sum(real(G_jk(OM_vect, x_j(s), x_k(s))- G_r_num(OM_vect, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2 + ...
            imag(G_jk(OM_vect, x_j(s), x_k(s))- G_r_num(OM_vect, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect, x_j(s), x_k(s)))); % find the index for which I have the maximum value of the FRF
om_i_0 = OM_vect(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(s), x_k(s)));
value2 = angle(G_jk(om_i_0+0.005, x_j(s), x_k(s)));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m);

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = imag(max(G_jk(OM_vect, x_j(s), x_k(s))))*1i * 2*1i*csi_0*om_i_0^2;

% Residuals = 0
% Vector of initial guesses
x0 = [om_i_0; csi_0; A_jk_i_0; 0; 0; 0; 0];


X = lsqnonlin(err,x0, [],[],[]);

figure (5+s)
subplot(2,1,1)
hold on
semilogy(OM_vect/2/pi, abs(G_r_num(OM_vect, X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize',8)
subplot(2,1,2)
hold on
grid on
plot(OM_vect/2/pi, angle(G_r_num(OM_vect,X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize', 8)

end

% MODE 2 

x_k=[1.2;0.2;0.605;1];
x_j=[0.2;1.2;0.8;0.605];

for s = 1:length(x_j)

G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;

OM_vect_2 = 2*pi*(27.8:0.005:28.8);

err = @ (X) sum(real(G_jk(OM_vect_2, x_j(s), x_k(s))- G_r_num(OM_vect_2, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2 + ...
            imag(G_jk(OM_vect_2, x_j(s), x_k(s))- G_r_num(OM_vect_2, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect_2, x_j(s), x_k(s)))); % find the index for which I have the maximum value of the FRF
om_i_0 = OM_vect_2(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(s), x_k(s)));
value2 = angle(G_jk(om_i_0+0.005, x_j(s), x_k(s)));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m);

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = imag(max(G_jk(OM_vect_2, x_j(s), x_k(s))))*1i * 2*1i*csi_0*om_i_0^2;

% Residuals = 0
% Vector of initial guesses
x0 = [om_i_0; csi_0; A_jk_i_0; 0; 0; 0; 0];


X = lsqnonlin(err,x0, [],[],[]);

figure (5+s)
subplot(2,1,1)
hold on
semilogy(OM_vect_2/2/pi, abs(G_r_num(OM_vect_2, X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize',8)
subplot(2,1,2)
hold on
grid on
plot(OM_vect_2/2/pi, angle(G_r_num(OM_vect_2,X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize', 8)

end

% MODE 3

x_k=[1.2;0.2;0.605;1];
x_j=[0.2;1.2;0.8;0.605];

for s = 1:length(x_j)

G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;

OM_vect_3 = 2*pi*(78.5:0.005:79.5);

err = @ (X) sum(real(G_jk(OM_vect_3, x_j(s), x_k(s))- G_r_num(OM_vect_3, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2 + ...
            imag(G_jk(OM_vect_3, x_j(s), x_k(s))- G_r_num(OM_vect_3, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect_3, x_j(s), x_k(s)))); % find the index for which I have the maximum value of the FRF
om_i_0 = OM_vect_3(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(s), x_k(s)));
value2 = angle(G_jk(om_i_0+0.005, x_j(s), x_k(s)));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m);

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = imag(max(G_jk(OM_vect_3, x_j(s), x_k(s))))*1i * 2*1i*csi_0*om_i_0^2;

% Residuals = 0
% Vector of initial guesses
x0 = [om_i_0; csi_0; A_jk_i_0; 0; 0; 0; 0];


X = lsqnonlin(err,x0, [],[],[]);

figure (5+s)
subplot(2,1,1)
hold on
semilogy(OM_vect_3/2/pi, abs(G_r_num(OM_vect_3, X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize',8)
subplot(2,1,2)
hold on
grid on
plot(OM_vect_3/2/pi, angle(G_r_num(OM_vect_3,X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize', 8)

end

% MODE 4

x_k=[1.2;0.2;0.605;1];
x_j=[0.2;1.2;0.8;0.605];

for s = 1:length(x_j)

G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;

OM_vect_4 = 2*pi*(154.03:0.005:155.5);

err = @ (X) sum(real(G_jk(OM_vect_4, x_j(s), x_k(s))- G_r_num(OM_vect_4, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2 + ...
            imag(G_jk(OM_vect_4, x_j(s), x_k(s))- G_r_num(OM_vect_4, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect_4, x_j(s), x_k(s)))); % find the index for which I have the maximum value of the FRF
om_i_0 = OM_vect_4(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(s), x_k(s)));
value2 = angle(G_jk(om_i_0+0.005, x_j(s), x_k(s)));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m);

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = imag(max(G_jk(OM_vect_4, x_j(s), x_k(s))))*1i * 2*1i*csi_0*om_i_0^2;

% Residuals = 0
% Vector of initial guesses
x0 = [om_i_0; csi_0; A_jk_i_0; 0; 0; 0; 0];


X = lsqnonlin(err,x0, [],[],[]);

figure (5+s)
subplot(2,1,1)
hold on
semilogy(OM_vect_4/2/pi, abs(G_r_num(OM_vect_4, X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize',8)
subplot(2,1,2)
hold on
grid on
plot(OM_vect_4/2/pi, angle(G_r_num(OM_vect_4,X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize', 8)

end









