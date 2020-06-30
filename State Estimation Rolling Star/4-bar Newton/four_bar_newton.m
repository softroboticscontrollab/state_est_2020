%% RISS 2020 4-bar Newton %%

clc
close all
clear variables

%% enter knowns
K1 = 1000;
K2 = 1000;
m1 = 7;
m2 = 7;
m3 = 7;
m1p = 7;
g = 9.81;
a1 = 1;
a2 = 1;
a3 = 1;
a1p = 1;
theta0=pi/2;
param = [K1,K2,m1,m2,m3,m1p,g,a1,a2,a3,a1p,theta0];

%% Newton Raphson
tol=0.001;
z1=pi/3;
z2=pi/3;
z3=pi/3;
z4=pi/3;
z5=pi/3;
zi=[z1;z2;z3;z4;z5];
maxIter=100;
toggle=1;
fprintf('Case 1 Iteration:\n');
[soln,er_est1]=func_MDnewton(@resid_vec_four,@dRdx_four,zi,tol,maxIter,toggle,param);

% x(1) = 0;
% y(1) = 0;
% x(2) = a1*cos(soln(1));
% y(2) = a1*sin(soln(1));
% x(3) = x(2)+a2*cos(soln(1)+soln(2));
% y(3) = y(2)+a2*sin(soln(1)+soln(2));
% soln_d = soln * 180/pi
% 
% plot(x,y,'k-')
% size = 2;
% axis([-size size -size size])
% axis square