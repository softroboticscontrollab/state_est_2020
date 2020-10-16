%% RISS 2020 statics 4-bar numeric solver %%
% Using statics, determine the necessary angles to fully define the shape
% of the robot.  This script utilizes the equler lagrange approach

clc
clear all
close all 

%% enter parameters
a1=2;
a2=3;
a3=2;
a4=4;
m1=1;
m2=1;
m3=1;
g=9.81;
k2=1000;
k3=1000;
t0=pi/2;

param = [a1 a2 a3 a4 m1 m2 m3 g k2 k3 t0]; %create a vector of all paramters

%% newton raphson 
%inputs and functions
R=@resid_single_pend; %residual function
dRdx=@dresid_single_pend; %slope of residual
t1i=pi/2; %initial guess for solution
tol=0.0001; %solution tolerance
maxIter=100; %max iteratios to find solution
toggle=1; %1 prints guesses

%Run numeric solver
[t1_sol,er_est]=func_newton(R,dRdx,t1i,tol,maxIter,toggle,param);

%% plot results
%% Create corner vector
%find the numeric value of t2 and t3 
x = sqrt(a1^2 + a4^2 - 2*a1*a4*cos(pi-t1));
alpha = asin((a4*sin(t1))/x);
beta = acos((x^2 + a2^2 - a3^2)/(2*x*a2));
zeta = asin((x*sin(beta))/a3);
t2_sol = simplify(pi - alpha - beta);
t3_sol = simplify(pi-zeta);

x(1) = 0;
y(1) = 0;
x(2) = a1*cos(t1_sol);
y(2) = a1*sin(t1_sol);
x(3) = x(2) + a2*cos(t2);
y(3) = y(2) + a2*sin(t2);
x(4) = -a4;
y(4) = 0);
x(5) = 0;
y(5) = 0;

plot(x,y)
axis equal
axis([-7 7 -7 7])



