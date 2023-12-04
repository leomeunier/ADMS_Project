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

% Searching for natural frequencies, so when determinant([H]) = 0

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

% Solving H*z = 0
% Since I found the first 4 omega0 that make det([H]) = 0, now the rows of
% [H] are not linearly independent anymore. So I remove the first row and I
% obtain H_new. While doing that I removed also the first column.
% So now I consider it and call it N.
% I obtained a new system to solve: H_new * z = - N
% Note that the first row of z is still full of ones because is the row I
% considered as known when at the beginning I've removed the first row of
% the matrix [H].

for i = 1:length(omega0)
    Hgamma = H(gamma0(i));
    H_new = Hgamma(2:end, 2:end);
    N = Hgamma(2:end,1);
    k = H_new\(-N);
    z(2:end, i) = k; % I found the coefficients needed to evaluate the transverse vibration w.
    % z is a matrix in which each column is a vector with the coefficients
    % for each natural frequency.
end

% I have w = phi(x) * G(t), but I only want the mode shapes so phi(x).

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


%% Numerical Frequency Response Function - Comparison with the one on slide 12

% MODE 1 

x_j = 0.2;
x_k = 1.2; 

G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;

OM_vect = 2*pi*(4:0.005:5); 

err = @ (X) sum(real(G_jk(OM_vect, x_j, x_k)- G_r_num(OM_vect, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2 + ...
            imag(G_jk(OM_vect, x_j, x_k)- G_r_num(OM_vect, X(1), X(2), X(3), X(4), X(5), X(6), X(7))).^2);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect, x_j, x_k))); % find the index for which I have the maximum value of the FRF
om_i_0 = OM_vect(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j, x_k));
value2 = angle(G_jk(om_i_0+0.005, x_j, x_k));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m);

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = imag(max(G_jk(OM_vect, x_j, x_k)))*1i * 2*1i*csi_0*om_i_0^2;

% Residuals = 0
% Vector of initial guesses
x0 = [om_i_0; csi_0; A_jk_i_0; 0; 0; 0; 0];


X = lsqnonlin(err,x0, [],[],[]);

figure (6)
subplot(2,1,1)
hold on
semilogy(OM_vect/2/pi, abs(G_r_num(OM_vect, X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize',8)
subplot(2,1,2)
hold on
grid on
plot(OM_vect/2/pi, angle(G_r_num(OM_vect,X(1), X(2), X(3), X(4), X(5), X(6), X(7))),'ro','MarkerSize', 8)


%% FREQUENCY RESPONSE FUNCTIONS for a given set of experimental FRFs

% First mode
% 3 FRFs, to compare with slide 13

x_j = [0.2, 0.42, 0.9];
x_k = 1.2;

% G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             %(Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;
G_r_num1 = @ (OM, om_i, csi_i, xx) (xx(:,1)./(-OM.^2+1i*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (xx(:,2)+1i*xx(:,4))./OM.^2 + xx(:,3)+1i*xx(:,5));
% By doing this I obtain:
% xx = [A_jk_1, Rlow_re1, Rlow_im1, Rhigh_re1, Rlow_im1]
%      ....
%      [A_jk_r, Rlow_rer, Rlow_imr, Rhigh_rer, Rlow_imr]
% r = number of x_j x_k couples

OM_vect = 2*pi*(4:0.005:5);
l = length(x_j);
err_r = @(x)sum(sum(real(G_jk(OM_vect,x_j,x_k) - G_r_num1(OM_vect, x(1,1),x(1,2), x(:,3:7))).^2 + ...
    imag(G_jk(OM_vect,x_j,x_k) - G_r_num1(OM_vect,x(1,1),x(1,2), x(:,3:7))).^2));

x_0 = zeros(length(x_j), 7);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect, x_j(1), x_k))); % find the index for which I have the maximum value of the FRF
% Since omegai is the same for each FRF, I just consider 1.
om_i_0 = OM_vect(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(1), x_k));
value2 = angle(G_jk(om_i_0+0.005, x_j(1), x_k));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m); % Just one because is the same for all FRFs

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = zeros(length(x_j), 1);

for ii = 1:length(x_j)
    A_jk_i_0(ii) = imag(max(G_jk(OM_vect, x_j(ii),x_k)))*1i * 2*1i*csi_0*om_i_0^2;
end

% Residuals = 0

% Vector of initial guesses
x_0(:,1) = om_i_0;
x_0(:,2) = csi_0;
x_0(:,3) = A_jk_i_0;

X = lsqnonlin(err_r,x_0,[],[],[]);

% I compute another graph for G_jk in order to compare this set of FRFs to
% the one found in the previous step
OM_vect2 = 2*pi*(0:0.02:200);

figure (7)
subplot(2,1,1)
semilogy(OM_vect2/2/pi, abs(G_jk(OM_vect2,x_j1,x_k1)), 'LineWidth', 2)
ylim([1e-6 0.3e-1])
xlabel('Frequency [Hz]')
ylabel('[m/N]'); 
legend('G_{jk} ( \Omega ) = x_{j}/F_{k}, x_{j} = 0.2 m , x_{k} = 1.2 m')
grid on; hold on;

subplot(2,1,2)
plot(OM_vect2/2/pi, angle(G_jk(OM_vect2,x_j1,x_k1)), 'LineWidth',2)
grid on

% Now I plot my set of FRFS

figure (7)
subplot(2,1,1)
hold on
semilogy(OM_vect/2/pi, abs(G_r_num1(OM_vect, X(1,1),X(1,2), X(1,3:7))), 'go','MarkerSize',8) %I take only the first row because I'm comparing it with
% G_jk which had 0.2 as x_j

%Plot of initial guesses
%semilogy(OM_vect/2/pi, abs(G_r_num1(OM_vect, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

subplot(2,1,2)
hold on
grid on
plot(OM_vect/2/pi, angle(G_r_num1(OM_vect, X(1,1),X(1,2), X(1,3:7))),'go','MarkerSize', 8)

% Plot of initial guesses
%plot(OM_vect/2/pi, angle(G_r_num1(OM_vect, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

% Plot the identified mode shapes 
M = [];
for ii = 1:length(x_j)
    [~, I] = max(abs(G_r_num1(OM_vect, X(1,1),X(1,2), X(ii,3:7))));
    M = [M; imag(G_r_num1(OM_vect(I), X(1,1),X(1,2), X(ii,3:7)))];
end

M = M * phi1(0.2)/M(1)/(max(abs(phi1(distance))));
figure (2)
plot(x_j,M,'ro', 'MarkerSize', 8, 'LineWidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second mode
% 3 FRFs, to compare with slide 13

x_j = [0.2, 0.42, 0.9];
x_k = 1.2;

% G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             %(Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;
G_r_num1 = @ (OM, om_i, csi_i, xx) (xx(:,1)./(-OM.^2+1i*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (xx(:,2)+1i*xx(:,4))./OM.^2 + xx(:,3)+1i*xx(:,5));
% By doing this I obtain:
% xx = [A_jk_1, Rlow_re1, Rlow_im1, Rhigh_re1, Rlow_im1]
%      ....
%      [A_jk_r, Rlow_rer, Rlow_imr, Rhigh_rer, Rlow_imr]
% r = number of x_j x_k couples

OM_vect_2 = 2*pi*(27.8:0.005:28.7);
l = length(x_j);
err_r = @(x)sum(sum(real(G_jk(OM_vect_2,x_j,x_k) - G_r_num1(OM_vect_2, x(1,1),x(1,2), x(:,3:7))).^2 + ...
    imag(G_jk(OM_vect_2,x_j,x_k) - G_r_num1(OM_vect_2,x(1,1),x(1,2), x(:,3:7))).^2));

x_0 = zeros(length(x_j), 7);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect_2, x_j(1), x_k))); % find the index for which I have the maximum value of the FRF
% Since omegai is the same for each FRF, I just consider 1.
om_i_0 = OM_vect_2(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(1), x_k));
value2 = angle(G_jk(om_i_0+0.005, x_j(1), x_k));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m); % Just one because is the same for all FRFs

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = zeros(length(x_j), 1);

for ii = 1:length(x_j)
    A_jk_i_0(ii) = imag(max(G_jk(OM_vect_2, x_j(ii),x_k)))*1i * 2*1i*csi_0*om_i_0^2;
end

% Residuals = 0

% Vector of initial guesses
x_0(:,1) = om_i_0;
x_0(:,2) = csi_0;
x_0(:,3) = A_jk_i_0;

X = lsqnonlin(err_r,x_0,[],[],[]);

figure (7)
subplot(2,1,1)
hold on
semilogy(OM_vect_2/2/pi, abs(G_r_num1(OM_vect_2, X(1,1),X(1,2), X(1,3:7))), 'go','MarkerSize',8) %I take only the first row because I'm comparing it with
% G_jk which had 0.2 as x_j

%Plot of initial guesses
%semilogy(OM_vect/2/pi, abs(G_r_num1(OM_vect_2, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

subplot(2,1,2)
hold on
grid on
plot(OM_vect_2/2/pi, angle(G_r_num1(OM_vect_2, X(1,1),X(1,2), X(1,3:7))),'go','MarkerSize', 8)

% Plot of initial guesses
% plot(OM_vect/2/pi, angle(G_r_num1(OM_vect, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

% Plot the identified mode shapes 
M = [];
for ii = 1:length(x_j)
    [~, I] = max(abs(G_r_num1(OM_vect_2, X(1,1),X(1,2), X(ii,3:7))));
    M = [M; imag(G_r_num1(OM_vect_2(I), X(1,1),X(1,2), X(ii,3:7)))];
end

M = M*phi2(0.2)/M(1)/(max(abs(phi2(distance))));
figure (3)
plot(x_j,M,'ro', 'MarkerSize', 8, 'LineWidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Third mode
% 3 FRFs, to compare with slide 13

x_j = [0.2, 0.42, 0.9];
x_k = 1.2;

% G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             %(Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;
G_r_num1 = @ (OM, om_i, csi_i, xx) (xx(:,1)./(-OM.^2+1i*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (xx(:,2)+1i*xx(:,4))./OM.^2 + xx(:,3)+1i*xx(:,5));
% By doing this I obtain:
% xx = [A_jk_1, Rlow_re1, Rlow_im1, Rhigh_re1, Rlow_im1]
%      ....
%      [A_jk_r, Rlow_rer, Rlow_imr, Rhigh_rer, Rlow_imr]
% r = number of x_j x_k couples

OM_vect_3 = 2*pi*(78.5:0.005:79.5);
l = length(x_j);
err_r = @(x)sum(sum(real(G_jk(OM_vect_3,x_j,x_k) - G_r_num1(OM_vect_3, x(1,1),x(1,2), x(:,3:7))).^2 + ...
    imag(G_jk(OM_vect_3,x_j,x_k) - G_r_num1(OM_vect_3,x(1,1),x(1,2), x(:,3:7))).^2));

x_0 = zeros(length(x_j), 7);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect_3, x_j(1), x_k))); % find the index for which I have the maximum value of the FRF
% Since omegai is the same for each FRF, I just consider 1.
om_i_0 = OM_vect_3(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(1), x_k));
value2 = angle(G_jk(om_i_0+0.005, x_j(1), x_k));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m); % Just one because is the same for all FRFs

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = zeros(length(x_j), 1);

for ii = 1:length(x_j)
    A_jk_i_0(ii) = imag(max(G_jk(OM_vect_3, x_j(ii),x_k)))*1i * 2*1i*csi_0*om_i_0^2;
end

% Residuals = 0

% Vector of initial guesses
x_0(:,1) = om_i_0;
x_0(:,2) = csi_0;
x_0(:,3) = A_jk_i_0;

X = lsqnonlin(err_r,x_0,[],[],[]);

figure (7)
subplot(2,1,1)
hold on
semilogy(OM_vect_3/2/pi, abs(G_r_num1(OM_vect_3, X(1,1),X(1,2), X(1,3:7))), 'go','MarkerSize',8) %I take only the first row because I'm comparing it with
% G_jk which had 0.2 as x_j

%Plot of initial guesses
%semilogy(OM_vect/2/pi, abs(G_r_num1(OM_vect_2, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

subplot(2,1,2)
hold on
grid on
plot(OM_vect_3/2/pi, angle(G_r_num1(OM_vect_3, X(1,1),X(1,2), X(1,3:7))),'go','MarkerSize', 8)

% Plot of initial guesses
% plot(OM_vect/2/pi, angle(G_r_num1(OM_vect, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

% Plot the identified mode shapes 
M = [];
for ii = 1:length(x_j)
    [~, I] = max(abs(G_r_num1(OM_vect_3, X(1,1),X(1,2), X(ii,3:7))));
    M = [M; imag(G_r_num1(OM_vect_3(I), X(1,1),X(1,2), X(ii,3:7)))];
end

M = M*phi3(0.2)/M(1)/(max(abs(phi3(distance))));
figure (4)
plot(x_j,M,'ro', 'MarkerSize', 8, 'LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fourth mode
% 3 FRFs, to compare with slide 13

x_j = [0.2, 0.42, 0.9];
x_k = 1.2;

% G_r_num =  @ (OM, om_i, csi_i, A_jk_i, Rlow_re, Rlow_im, Rhigh_re, Rhigh_im)  A_jk_i./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + ...
             %(Rlow_re+1j*Rlow_im)./OM.^2 + Rhigh_re+1j*Rhigh_im;
G_r_num1 = @ (OM, om_i, csi_i, xx) (xx(:,1)./(-OM.^2+1i*2*csi_i.*om_i.*OM+om_i.^2) + ...
             (xx(:,2)+1i*xx(:,4))./OM.^2 + xx(:,3)+1i*xx(:,5));
% By doing this I obtain:
% xx = [A_jk_1, Rlow_re1, Rlow_im1, Rhigh_re1, Rlow_im1]
%      ....
%      [A_jk_r, Rlow_rer, Rlow_imr, Rhigh_rer, Rlow_imr]
% r = number of x_j x_k couples

OM_vect_4 = 2*pi*(154.03:0.005:155.5);
l = length(x_j);
err_r = @(x)sum(sum(real(G_jk(OM_vect_4,x_j,x_k) - G_r_num1(OM_vect_4, x(1,1),x(1,2), x(:,3:7))).^2 + ...
    imag(G_jk(OM_vect_4,x_j,x_k) - G_r_num1(OM_vect_4,x(1,1),x(1,2), x(:,3:7))).^2));

x_0 = zeros(length(x_j), 7);

% Initial guesses

% omegai
[~,om_index] = max(abs(G_jk(OM_vect_4, x_j(1), x_k))); % find the index for which I have the maximum value of the FRF
% Since omegai is the same for each FRF, I just consider 1.
om_i_0 = OM_vect_4(om_index);

% Damping ratio
% csi = - 1 / ( om_i_0 * dphi/dbigomega|om_i_0 );
value1 = angle(G_jk(om_i_0-0.005, x_j(1), x_k));
value2 = angle(G_jk(om_i_0+0.005, x_j(1), x_k));
m = (value2-value1)/((om_i_0 + 0.005)-(om_i_0-0.005));
csi_0 = -1/(om_i_0*m); % Just one because is the same for all FRFs

% A_jk
% At the peak we can approximate G_num as Ajk/(2j*csi0*om_i_0^2)
A_jk_i_0 = zeros(length(x_j), 1);

for ii = 1:length(x_j)
    A_jk_i_0(ii) = imag(max(G_jk(OM_vect_4, x_j(ii),x_k)))*1i * 2*1i*csi_0*om_i_0^2;
end

% Residuals = 0

% Vector of initial guesses
x_0(:,1) = om_i_0;
x_0(:,2) = csi_0;
x_0(:,3) = A_jk_i_0;

X = lsqnonlin(err_r,x_0,[],[],[]);

figure (7)
subplot(2,1,1)
hold on
semilogy(OM_vect_4/2/pi, abs(G_r_num1(OM_vect_4, X(1,1),X(1,2), X(1,3:7))), 'go','MarkerSize',8) %I take only the first row because I'm comparing it with
% G_jk which had 0.2 as x_j

%Plot of initial guesses
%semilogy(OM_vect/2/pi, abs(G_r_num1(OM_vect_2, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

subplot(2,1,2)
hold on
grid on
plot(OM_vect_4/2/pi, angle(G_r_num1(OM_vect_4, X(1,1),X(1,2), X(1,3:7))),'go','MarkerSize', 8)

% Plot of initial guesses
% plot(OM_vect/2/pi, angle(G_r_num1(OM_vect, x_0(1,1),x_0(1,2), x_0(1,3:7))),'bx', 'MarkerSize', 8)

% Plot the identified mode shapes 
M = [];
for ii = 1:length(x_j)
    [~, I] = max(abs(G_r_num1(OM_vect_4, X(1,1),X(1,2), X(ii,3:7))));
    M = [M; imag(G_r_num1(OM_vect_4(I), X(1,1),X(1,2), X(ii,3:7)))];
end

M = M*phi4(0.2)/M(1)/(max(abs(phi4(distance))));
figure (5)
plot(x_j,M,'ro', 'MarkerSize', 8, 'LineWidth',2)



%% NOTES

% Consider more than 3 FRFs and maybe change position of the force
% Fix the legends 
% Table comparing model with identified as shown in slide 13









