
close all

deg2rad = pi/180;
alp = 120*deg2rad;
bet = 120*deg2rad;
del = 2*pi - alp - bet;

gam = 100*deg2rad;

a = 1; 
b = 0.5;    c = 0.5;
d = a;
e = c;



D = [0,0];
E = D + e*[cos(del), sin(del)];
C = D + [d,0];
B = C + c*[cos(pi-gam), sin(pi-gam)];
A = B - b*[sin(alp-gam+del-pi/2), -cos(alp-gam+del-pi/2)];
F = A - a*[cos(gam-del), -sin(gam-del)]




% P6 = P5 + f*[cos(E-pi+D), sin(E-pi+D)];
% ;
% P22 = P1 + b*[sin(A-C+D-pi/2), -cos(A-C+D-pi/2)];


PP = [B;C;D;E;F;A;B];

plot(PP(:,1), PP(:,2), 'o-')

abs(sin(alp)*sin(eps)/sin(bet)/sin(zet))