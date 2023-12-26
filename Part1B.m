%% Part B

clc
close all
clear all

%% Load data
[filename, pathname] = uigetfile;
cd(pathname); 
load(filename);

%% Experimental FRF
figure (1)
subplot(3,1,1)
for ii = 1:size(frf,2)
    semilogy(freq, abs(frf(:,ii)),'LineWidth', 1.25)
    grid on
    hold on
end
xlabel('Frequency [Hz]')
ylabel('Amplitude [m/s^2/N]')
title('FrF amplitude')
legend('1','2','3','4','5','6','7','8','9','10','11','12')
subplot(3,1,2)
for ii = 1:size(frf,2)
    plot(freq, angle(frf(:,ii)),'LineWidth', 1.25)
    grid on
    hold on
end
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
title('FrF phase')

subplot(3,1,3)
for ii = 1:size(frf,2)
    semilogy(freq, cohe(:,ii),'LineWidth', 1.25)
    grid on
    hold on
end
xlabel('Frequency [Hz]')
ylabel('Amplitude [m/s^2/N]')
title('FrF coherence')

frf1 = frf(:,1);

maxf = 0;
maxfrf1 = 0;

for ii = 2:length(frf1)
    if abs(frf1(ii)) > abs(frf1(ii-1)) && abs(frf1(ii)) > abs(frf1(ii+1))
        maxfrf1 = abs(frf1(ii));
        maxf = freq(ii);
    end

    if maxf > 666 && maxf < 669
        nat_freq1 = maxf;
    end

    if maxf > 1615 && maxf < 1640
        nat_freq2 = maxf;
    end
end

freq0 = [nat_freq1 nat_freq2];
omega0 = freq0*2*pi;



%% Modal parameters identification

%% Estimations for first peak


freq1 = 6.499913334488874e+02;
freq2 = 6.796576045652725e+02;

for ii = 1:length(freq)
    if freq(ii) == freq1
        index1 = ii;
    end

    if freq(ii) == freq2
        index2 = ii;
    end
end

OM_vect = 2*pi*(freq(index1):(freq(index1)-freq(index1-1)):freq(index2));
min_frf = frf(index1:index2,:);

w1= 661.97*2*pi;
w2= 672.71*2*pi;
csi_0= (w2^2-w1^2)/(4*omega0(1)^2);

% A
A0 = zeros(1,12);
for ii = 1:12
    A0(ii) = imag(max(frf((index1:index2),ii))).*1i.*2.*1i.*omega0(1).^2.*csi_0;
end

%Vector of initial guesses
x_0 = zeros(12,7);
x_0(:,1) = omega0(1);
x_0(:,2) = csi_0;
x_0(:,3) = A0;

G_r_num1 =  @(OM, om_i,csi_i, xx)conj(xx(:,1)./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + (xx(:,2)+1j*xx(:,4))./OM.^2 + xx(:,3)+1j*xx(:,5))';

err_r = @(x)sum(sum(real(min_frf - G_r_num1(OM_vect, x(1,1),x(1,2), x(:,3:7))).^2 + imag(min_frf- G_r_num1(OM_vect,x(1,1),x(1,2), x(:,3:7))).^2));

X = lsqnonlin(err_r,x_0,[],[],[]);

for ii = 1:12
    figure(ii+1)
    subplot(2,1,1)
    semilogy(freq, abs(frf(:,ii)),'LineWidth', 1.25)
    hold on
    grid on
    semilogy(OM_vect/2/pi, abs(G_r_num1(OM_vect, X(1,1),X(1,2), X(ii,3:7))),'go', 'MarkerSize',8)
    grid on
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [m/s^2/N]')
    title('FrF amplitude')
    subplot(2,1,2)
    plot(freq, angle(frf(:,ii)),'LineWidth', 1.25)
    hold on
    grid on
    plot(OM_vect/2/pi, angle(G_r_num1(OM_vect, X(1,1),X(1,2), X(ii,3:7))),'go','MarkerSize', 8)
    grid on
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Phase [rad]')
    title('FrF phase')
end

%% MODE SHAPE for first natural frequency

theta = [pi/2:15*pi/180:3*pi/2];
Peaks = [];
for ii = 1:12
    [~, I] = max(abs(G_r_num1(OM_vect, X(1,1),X(1,2), X(ii,3:7))));
    Peaks = [Peaks; imag(G_r_num1(OM_vect(I), X(1,1),X(1,2), X(ii,3:7)))];
end

Peaks = [Peaks; Peaks(1)];
figure(14)
polarplot(theta,Peaks+5)
hold on
polarplot(theta,-(Peaks+5))
hold on
theta2 = [pi/2:15*pi/180:450*pi/180];
polarplot(theta2,5*ones(25,1))
legend('Tire axial displacement', 'Symmetric', 'Undeformed tire')
 










%% Estimations for second peak

index1 = 4852;
index2 = 4900;

OM_vect = 2*pi*(freq(index1):(freq(index1)-freq(index1-1)):freq(index2));
min_frf = frf(index1:index2-1,:);


w1= 1615.24*2*pi;
w2= 1635.09*2*pi;
damp2= (w2^2-w1^2)/(4*omega0(2)^2);

% A
A0 = zeros(1,12);
for ii = 1:12
    A0(ii) = imag(max(frf((index1:index2),ii))).*1i.*2.*1i.*omega0(2).^2.*damp2;
end

%Vector of initial guesses
x_0 = zeros(12,7);
x_0(:,1) = omega0(2);
x_0(:,2) = damp2;
x_0(:,3) = A0;

G_r_num1 =  @(OM, om_i,csi_i, xx)conj(xx(:,1)./(-OM.^2+1j*2*csi_i.*om_i.*OM+om_i.^2) + (xx(:,2)+1j*xx(:,4))./OM.^2 + xx(:,3)+1j*xx(:,5))';

err_r = @(x)sum(sum(real(min_frf - G_r_num1(OM_vect, x(1,1),x(1,2), x(:,3:7))).^2 + imag(min_frf- G_r_num1(OM_vect,x(1,1),x(1,2), x(:,3:7))).^2));

X = lsqnonlin(err_r,x_0,[],[],[]);

for ii = 1:12
    figure(ii+1)
    subplot(2,1,1)
    hold on
    grid on
    semilogy(OM_vect/2/pi, abs(G_r_num1(OM_vect, X(1,1),X(1,2), X(ii,3:7))),'go', 'MarkerSize',8)
    grid on
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [m/s^2/N]')
    title('FrF amplitude')
    subplot(2,1,2)
    hold on
    grid on
    plot(OM_vect/2/pi, angle(G_r_num1(OM_vect, X(1,1),X(1,2), X(ii,3:7))),'go','MarkerSize', 8)
    grid on
    hold on
    xlabel('Frequency [Hz]')
    ylabel('Phase [rad]')
    title('FrF phase')
end

%% MODE SHAPE for second natural frequency

theta = [pi/2:15*pi/180:3*pi/2];
Peaks = [];
for ii = 1:12
    [~, I] = max(abs(G_r_num1(OM_vect, X(1,1),X(1,2), X(ii,3:7))));
    Peaks = [Peaks; imag(G_r_num1(OM_vect(I), X(1,1),X(1,2), X(ii,3:7)))];
end

Peaks = [Peaks; Peaks(1)];
figure(15)
polarplot(theta,Peaks+5)
hold on
%polarplot(theta,-(Peaks+5))
hold on
theta2 = [pi/2:15*pi/180:450*pi/180];
polarplot(theta2,5*ones(25,1))
legend('Tire axial displacement', 'Undeformed tire')