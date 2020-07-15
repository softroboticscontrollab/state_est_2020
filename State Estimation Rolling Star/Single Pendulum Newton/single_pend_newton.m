%% RISS 2020 Single Pendulum %%
% given the spring properties, mass, geometry, and angle of a single
% pendulum, calculate the angle at which the pendulum will settle.  This
% uses the newton raphson numeric solver

clc
clear all
close all 

%% enter parameters
K=100; %curvature (m^-1)
m=7; %mass (kg)
g=9.81; %acceleration due to gravity (m/s/s)
l=1; %bar length (m)
theta0=pi/4; %unstretched angle (rad)
param = [K,m,g,l,theta0]; %create a vector of all paramters

%newton raphson inputs and functions
R=@resid_single_pend; %residual function
dRdx=@dresid_single_pend; %slope of residual
xi=0; %initial guess for solution
tol=0.0001; %solution tolerance
maxIter=100; %max iteratios to find solution
toggle=1; %1 prints guesses

%% Run numeric solver
[soln,er_est]=func_newton(R,dRdx,xi,tol,maxIter,toggle,param);

%% Plot results
x(1) = 0;
y(1) = 0;
x(2) = l*cos(soln);
y(2) = l*sin(soln);
soln_d = soln * 180/pi

p1 = plot(x,y,'k-',x(end),y(end),'rx')
xlabel('x-pos')
ylabel('y-pos')
linkLab = sprintf('link (angle = %.1f)',soln_d)
legend(linkLab,'COM','Location','south')
width = 2;
set(p1(1),{'LineWidth'},{width})

%max axis square to allow for representative viewing
size = 1;
axis([-size size -size size])
axis square
grid on