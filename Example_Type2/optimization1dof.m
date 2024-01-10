
clc
close all
clearvars



global A B C D E F a b c d e f

deg2rad = pi/180;
A = 140*deg2rad;
B = 140*deg2rad;
C = 90*deg2rad;

[bopt, fopt] = fminsearch(@fun, [0.5, 0.5])

P4 = [0,0];
P5 = P4 + e*[cos(D), sin(D)];
P3 = P4 + [d,0];
P2= P3+ c*[cos(pi-C), sin(pi-C)];
P1= P2- b*[sin(A-C+D-pi/2), -cos(A-C+D-pi/2)];
P6= P1- a*[cos(C-D), -sin(C-D)];


PP = [P2;P3;P4;P5;P6;P1;P2];
plot(PP(:,1), PP(:,2), 'o-')

[D,C,B,A,F,E]*180/pi

function [err] = fun(x)
global A B C D E F a b c d e f

b=x(1); c=x(2);

deg2rad = pi/180;
D = 2*pi - A - B;


a = 1; 
d = a;
e = c;


lb = - b * cos(A);
le = - e * cos(D);
la = a * cos(A+B+C);
lc = - c * cos(C);
lf = lc + la + le - a - lb;

hb = b * sin(A);
he = e * sin(D);
ha = a * sin(A+B+C);
hc = c * sin(C);
hf = hb+he+ha-hc;

f = (lf^2 + hf^2)^0.5;
F = pi-atan2(hf,lf);
E = 4 * pi - A - B - C - D - F ;


err = abs(abs(sin(A)*sin(E)/sin(B)/sin(F))-1);

% plot([b], [err], 'o')
% hold on

end



