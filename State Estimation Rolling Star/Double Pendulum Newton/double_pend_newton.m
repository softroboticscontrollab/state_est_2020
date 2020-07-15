%% RISS 2020 Double Pendulum %%
% given the spring properties, mass, geometry, and angle of a double
% pendulum, calculate the angle at which the pendulums will settle.  This
% uses the multi-dimensionsal newton raphson numeric solver

clc
close all
clear variables

%% Enter knowns (arbitrary)
K1 = 3000; %curvature 1 (m^-1)
K2 = 25;  %curvature 2 (m^-1)
m1 = 7; %mass 1 (kg)
m2 = 7; %mass 2 (kg)
g = 9.81; %acceleration due to gravity (m/s/s)
L1 = 1; %length of bar 1
L2 = 1; %length of bar 2
theta0=pi/6; %unstretched angle of springs 
param = [K1,K2,m1,m2,g,L1,L2,theta0]; %create a vector of parameters

%% Newton Raphson
tol=0.001; %solution tolerance
x1=0; %initial guess for the first angle
x2=pi/2; %initial guess for the second angle
xi=[x1;x2]; %create a vector of guesses
maxIter=5000; %max iterations to find solution
toggle=1; %1 prints solution
fprintf('Case 1 Iteration:\n');
[soln,er_est1]=func_MDnewton(@resid_vec,@dRdx,xi,tol,maxIter,toggle,param);

%kineamtics for plot/COM calcs
x(1) = 0;
y(1) = 0;
x(2) = L1*cos(soln(1));
y(2) = L1*sin(soln(1));
x(3) = x(2)+L2*cos(soln(1)+soln(2));
y(3) = y(2)+L2*sin(soln(1)+soln(2));
soln_d = soln * 180/pi

%determine COM
xcom(1) = (L1/2)*cos(soln(1));
ycom(1) = (L1/2)*sin(soln(1));
xcom(2) = x(2)+(L2/2)*cos(soln(1)+soln(2));
ycom(2) = y(2)+(L2/2)*sin(soln(1)+soln(2));

%plot 
p1 = plot(x(1:2),y(1:2),'k-',x(2:3),y(2:3),'b-',xcom,ycom,'rx')
size = 2;
axis square
grid on
axis([-size size -size size])
xlabel('x-pos')
ylabel('y-pos')
link1Lab = sprintf('link 1 (angle = %.1f)',soln_d(1))
link2Lab = sprintf('link 2 (angle = %.1f)',soln_d(2))
legend(link1Lab,link2Lab,'COM','Location','south')
width = 2;
set(p1(1),{'LineWidth'},{width})
set(p1(2),{'LineWidth'},{width})
set(p1(3),{'MarkerSize'},{12})


