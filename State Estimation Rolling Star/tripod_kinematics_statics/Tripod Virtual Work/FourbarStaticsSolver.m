%% FourbarStaticsSolver.m
%   Authors: Dr. Eric Constans (Rose-Hulman)
%            Sam Alvares (Rose-Hulman)
%   Date: 11/11/2020
%
%   Description: this program calls a numeric solver to determine the pose
%   of the fourbar robot at static equilibrium. It uses a virtual work
%   approach (VWork.m)and a Newton Raphson Solver (NewtonRaphsonVWork.m).
%   The four bar mechansim has two torsional springs at the non-ground
%   joints.

%% Prep the workspace
clc
clear all
close all

%% Parameters
% bar lengths
% a1 = 4;
% a2 = 4;
% a3 = 3;
% a4 = 2;
% [a,b,c,d] = sicilano_to_constans(a1,a2,a3,a4)

a = 1;         % crank length (m)
b = 1;         % coupler length (m)
c = 1;         % rocker length (m)
d = 1;         % length between ground pins (m)

% masses (kg)
m2 = 1;  m3 = 1; m4 = 1;

% acceleration due to gravity (m/s/s)
g = 9.81;

% spring constans (Nm/rad)
kB = 10; kC = 10;

%% Determine the configuration at the initial guess theta2_0
% inital guess for theta2 (rad)
theta2_0 =  pi/2;

% determine other angles at initial guess of theta2_0 (rad)
[theta3_0,theta4_0,h,delta] = fourbar(theta2_0,a,b,c,d);

% determine location of pins at initial guess of theta2_0 (m)
xA_0 = [0;0];
xB_0 = a*[cos(theta2_0); sin(theta2_0)];
xC_0 = xB_0 + b*[cos(theta3_0); sin(theta3_0)];
xD_0 = [d; 0];
xD2_0 = xC_0 - c*[cos(theta4_0); sin(theta4_0)];

% net deflection constants (rad), uncomment thetaB=-pi/2 and thetaC = pi/2
% if you want to have the springs be unstretched at 90 deg
thetaB = -pi/2;
thetaC = pi/2;
% thetaB = theta3_0 - theta2_0;
% thetaC = theta4_0 - theta3_0;

% plot pose at initial theta2 guess
figure; hold on; grid off
p1 = plot([xA_0(1) xB_0(1)],[xA_0(2) xB_0(2)],'k','LineWidth',1);
plot([xB_0(1) xC_0(1)],[xB_0(2) xC_0(2)],'k','LineWidth',1)
plot([xC_0(1) xD_0(1)],[xC_0(2) xD_0(2)],'k','LineWidth',1)
plot([xC_0(1) xD2_0(1)],[xC_0(2) xD2_0(2)],'k','LineWidth',1)

%% Run Newton Raphson solver
% create vector of parameters
param = [a b c d m2 m3 m4 g kB kC thetaB thetaC];
tol=0.0000001; % solution tolerance
maxIter=10000; % max iteratios to find solution
toggle=1; %1 prints guesses

% run numeric solver
[theta2,er_est]=NewtonRaphsonVWork(@VWork,theta2_0,tol,maxIter,toggle,param);

%% Determine the final configuration
% determine final angles
[theta3,theta4,g,h] = fourbar(theta2,a,b,c,d);

% determine location of pins at final configuration (m)
xA = [0;0];
xB = a*[cos(theta2); sin(theta2)];
xC = xB + b*[cos(theta3); sin(theta3)];
xD = [d; 0];
xD2 = xC - c*[cos(theta4); sin(theta4)];

% plot final configuration
p2 = plot([xA(1) xB(1)],[xA(2) xB(2)],'r','LineWidth',1);
plot([xB(1) xC(1)],[xB(2) xC(2)],'r','LineWidth',1)
plot([xC(1) xD(1)],[xC(2) xD(2)],'r','LineWidth',1)
plot([xC(1) xD2(1)],[xC(2) xD2(2)],'r','LineWidth',1)
yline(0)
axis equal
legend([p1,p2],'Configuration @ \theta_{2_0}','Final configuration')
xlabel('x-position (m)')
ylabel('y-position (m)')

% print results in command window
fprintf('\nSolution Angles:\n')
fprintf('theta2: %1.4f rad\ntheta3: %1.4f rad\ntheta4: %1.4f rad\n',theta2,theta3,theta4)

% Check energy simplified
yg2i = (a/2)*sin(theta2_0);
yg3i = a*sin(theta2_0)+(b/2)*sin(theta3_0);
yg4i = c/2*sin(theta4_0);

yg2f = (a/2)*sin(theta2);
yg3f = a*sin(theta2)+(b/2)*sin(theta3);
yg4f = c/2*sin(theta4);

thetaB = -pi/2; thetaC = pi/2;
phiBf = theta3 - theta2 - thetaB;
phiCf = theta4 - theta3 - thetaC;
phiBi = theta3_0 - theta2_0 - thetaB;
phiCi = theta4_0 - theta3_0 - thetaC;

m2 = 1;  m3 = 1; m4 = 1; g = 9.81;
Ei = g*(m2*yg2i+m3*yg3i+m4*yg4i)+.5*10*(phiBi^2+phiCi^2);
Ef = g*(m2*yg2f+m3*yg3f+m4*yg4f)+.5*10*(phiBf^2+phiCf^2);
fprintf('\nInitiail energy: %.9f\n', Ei)
fprintf('Final energy: %.9f\n', Ef)

