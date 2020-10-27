%% RISS 2020 statics 4-bar numeric solver Constans Kin%%
% Using the solution of four_bar_statics_v4, find the pose of the 4-bar
% mechanism at static equilibrium

%prep the workspace
clc
clear all
close all

%% enter parameters
n = 4; %number bars
a1 = 5;
L2 = 5;
L3 = 5;
L4 = 5;
kappa2 = .1;
kappa3 = .4;
kappa4 = .6;
m2=1;
m3=1;
m4=1;
g=9.81;
K3=100;
Kdel=100;
t03=pi/2;
t0del=pi/2;

% determine bar lengths
K = [kappa2; kappa3; kappa4];
L = [L2; L3; L4];
a(1) = a1;
for i = 1:n-1
    a(i+1) = barcalc(K(i),L(i));
end

param = [a1 L2 L3 L4 kappa2 kappa3 kappa4 m2 m3 m4 g K3 Kdel t03 t0del]; %create a vector of all paramters

%% newton raphson
%inputs and functions
R=@resid_four_bar_statics_constans; %residual function
dRdx=@dresid_four_bar_statics_constans; %slope of residual
t2i=pi/1.9; %initial guess for solution
tol=0.00001; %solution tolerance
maxIter=1000; %max iteratios to find solution
toggle=1; %1 prints guesses

%Run numeric solver
[t2_sol,er_est]=func_newton(R,dRdx,t2i,tol,maxIter,toggle,param);

%% plot results
%find the numeric value of t2 and t3
f_squared = a(1)^2+a(2)^2-2*a1*a(2)*cos(t2_sol)
del = acos((a(3)^2+a(4)^2-f_squared)/(2*a(3)*a(4)));
g_dist = a(3)-a(4)*cos(del);
h = a(4)*sin(del);
r = a(1)-a(2)*cos(t2_sol);
s = a(2)*sin(t2_sol);
t3_sol = atan((h*r-g_dist*s)/(g_dist*r+h*s));
t4_sol = del+t3_sol;

x1_vec = [0;0;0];
x2_vec = [0;0;0];
x3_vec = [a(2)*cos(t2_sol);
          a(2)*cos(t2_sol);
          0];
x4_vec = [a(1);0;0];

% Create corner vector
x(1) = 0;
y(1) = 0;
x(2) = a(2)*cos(t2_sol);
y(2) = a(2)*sin(t2_sol);
x(3) = a(1)+a(4)*cos(t4_sol);
y(3) = a(4)*sin(t4_sol);
x(4) = a(1);
y(4) = 0;
x(5) = 0;
y(5) = 0;

plot(x,y)
axis equal
axis([-2 7 -2 7])

%% Verify solution with energy
%determine the location of the COM
% for i = 1:n-1
%     phi = L(i)*K(i)*.5;
%     A =(2*sin(phi))/(L(i)*K(i)^2);
%     B = cos(phi)/K(i);
%     C = A-B;
%     if i==1
%         P0Gx(1) = a(1)/2*T0_n{1}(1,1)+C*T0_n{1}(1,2);
%         P0Gy(1) = a(1)/2*T0_n{1}(2,1)+C*T0_n{1}(2,2);
%     else
%         P0Gx(i) = T0_n{i-1}(1,4)+a(i)/2*T0_n{i}(1,1)+C*T0_n{i}(1,2);
%         P0Gy(i) = T0_n{i-1}(2,4)+a(i)/2*T0_n{i}(2,1)+C*T0_n{i}(2,2);
%     end
% end
% hold on
% plot(P0Gx,P0Gy,'rx')
% 
% % Determine initial PE with unstretched spring lengths
% Uinitial = g*(m1*L1/2+m2*L2+m3*L3/2)
% 
% Ufinal = g*(m1*P0Gy(1)+m2*P0Gy(2)+m3*P0Gy(3))+.5*k2*(pi-t2_sol-t02)^2+.5*(pi-t3_sol-t03)^2

% determine initial pose
xi(1) = 0;
yi(1) = 0;
xi(2) = L(2)*cos(pi/2);
yi(2) = L(2)*sin(pi/2);
xi(3) = L(3);
yi(3) = L(2);
xi(4) = L(3);
yi(4) = 0;
xi(5) = 0;
yi(5) = 0;

% determine inital COM pos
P0Gxi = [0;L(3)/2;L(3)];
P0Gyi = [L(2)/2;L(2);L(2)/2];

% Plot inital pose
hold on 
plot(xi,yi,'--k')

%Plot initial COM position
hold on 
plot(P0Gxi,P0Gyi,'kx')


