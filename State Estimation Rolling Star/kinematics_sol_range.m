%% Kineamtics solution ranges %%
% Determine the range of solutions that the kinematic solution works
% graphically

clc
clear all
close all

%% Enter Constants
a1 = 3;
a2 = 2;
a3 = 3;
a1p = 2;
% t1 = 0;
t1p = [0:0.01:2*pi];

for i = 1:length(t1p)
    t1(i) = 0;
%% Solve kineamtics
x = sqrt(a1^2 + a1p^2 - 2*a1*a1p*cos(t1p(i)-t1(i)));
alpha = asin((a1p*sin(t1p(i)-t1(i)))/x);
beta = acos((x^2 + a2^2 - a3^2)/(2*x*a2));
% s_zeta = (x*sin(beta))/a3;
% c_zeta = (-x^2+a3^2+a2^2)/(2*a3*a2);
% zeta = atan2(s_zeta,c_zeta);
t2 = pi - alpha - beta;
% t3_final = pi-zeta;

%% Create corner vector
x(1) = 0;
y(1) = 0;
x(2) = a1*cos(t1(i));
y(2) = a1*sin(t1(i));
x(3) = x(2) + a2*cos(t1_sol+t2_sol);
y(3) = y(2) + a2*sin(t1_sol+t2_sol);
x(4) = a1p*cos(t1p(i));
y(4) = a1p*sin(t1p(i));
x(5) = 0;
y(5) = 0;

plot(x,y)
axis equal
axis([-7 7 -7 7])

pause(.005)

end


