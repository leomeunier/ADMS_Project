clear all

%properties
red = [20e-3 1e-3];
blue = [32e-3 1.5e-3];
green = [50e-3 2.5e-3];

E = 7e10;
rho = 2700;

colors = [ red ; blue ; green];
m = zeros(1,3);
EA = zeros(1,3);
EJ = zeros(1,3);
i = 0;

for i = 1:3
    m(i)= rho*pi*(colors(i,1)^2 - (colors(i,1)-colors(i,2))^2)/4;
    EA(i)= E*pi*(colors(i,1)^2 - (colors(i,1)-colors(i,2))^2)/4;
    EJ(i)= E*(pi/2)*((colors(i,1)/2)^4 - ((colors(i,1)-colors(i,2))/2)^4);
end

c2 = min(EJ./m);
c = sqrt(c2);
fmax = 200;
Lmax = sqrt((1/(1.5*2*fmax))*c);

A = [0.7 0];
B = [-0.42 0];
C = [0 0];
D = [0.56 0.42];
E = [0.52 0.54];
F = [0.49 0.63];
G = [-0.08 0.48];
H = [-0.12 0.72];



