%% RISS 2020 Single Pendulum %%

clc
clear all
close all 

%% enter parameters
K=100;
m=7;
g=9.81;
l=1;
theta0=pi/4;
param = [K,m,g,l,theta0];

R=@resid_single_pend;
dRdx=@dresid_single_pend;
xi=0;
tol=0.0001;
maxIter=100;
toggle=1;

[soln,er_est]=func_newton(R,dRdx,xi,tol,maxIter,toggle,param);

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

size = 1;
axis([-size size -size size])
axis square
grid on