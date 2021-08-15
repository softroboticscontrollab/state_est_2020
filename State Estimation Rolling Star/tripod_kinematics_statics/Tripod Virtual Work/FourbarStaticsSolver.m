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

%% Set configurations
gravity = false;
neutral_right_angles = true;
force_to_right = true;
test_case = 2;

%% Parameters
if (test_case == 1)
    a1 = 4;
    a2 = 4;
    a3 = 5;
    a4 = 2;
elseif (test_case == 2)
    a1 = 4;
    a2 = 5;
    a3 = 3;
    a4 = 6;
elseif (test_case == 3)
    a1 = 1;
    a2 = 1;
    a3 = 1;
    a4 = 1;
end

[a,b,c,d] = sicilano_to_constans(a1,a2,a3,a4);

% a = crank length (m)
% b = coupler length (m)
% c = rocker length (m)
% d = length between ground pins (m)

% Masses (kg)
m2 = 1;  m3 = 1; m4 = 1;

% Acceleration due to gravity (m/s/s)
g = 9.81;

% Spring constans (Nm/rad)
kB = 50; kC = 50;

% Force to right
Fmax = 10;
Fstep = 0.0001;
F = 0:Fstep:Fmax;

% Forced due to gravity
Fg2 = -m2*g;
Fg3 = -m3*g;
Fg4 = -m4*g;

% net deflection constants (rad), uncomment thetaB=-pi/2 and thetaC = pi/2
% if you want to have the springs be unstretched at 90 deg
if (neutral_right_angles)
    thetaB = -pi/2;
    thetaC = pi/2;
else % calibrated values to be found after completing calibration...
    %     thetaC = thetaC_calibrated;
    %     thetaC = thetaC_calibrated;
end


% plot pose at initial theta2 guess
figure; hold on; grid off

%% Run Newton Raphson solver
% inital guess for theta2 (rad)
theta2_0 =  pi/1.9;
tol=0.0000001; % solution tolerance
maxIter=10000; % max iteratios to find solution
toggle=0; %1 prints guesses

for i = 1:length(F)
    % create vector of parameters
    param = [a b c d m2 m3 m4 g kB kC thetaB thetaC F(i) Fg2 Fg3 Fg4 gravity force_to_right];
    % run numeric solver
    [theta2(i)]=NewtonRaphsonVWork(@VWork,theta2_0,tol,maxIter,toggle,param);
    % determine all bar angles
    [theta3(i),theta4(i),g,h] = fourbar(theta2(i),a,b,c,d);
end

% determine location of pins at final configuration (m)
% x1(1) = 0;
% y1(1) = 0;
% x1(2) = a*cos(theta2(1));
% y1(2) = a*sin(theta2(1));
% x1(3) = x1(2) + b*cos(theta3(1));
% y1(3) = y1(2) + b*sin(theta3(1));
% x1(4) = d;
% y1(4) = 0;
% 
% xend(1) = 0;
% yend(1) = 0;
% xend(2) = a*cos(theta2(end));
% yend(2) = a*sin(theta2(end));
% xend(3) = xend(2) + b*cos(theta3(end));
% yend(3) = yend(2) + b*sin(theta3(end));
% xend(4) = d;
% yend(4) = 0;
% 
% figure(2)
% plot(x1,y1,'r')
% hold on
% plot(xend,yend,'b')
% axis equal
% hold off

%% Energy Verification
% Find work due to force, F
if (force_to_right)
    Wf_sum = 0;
    for i = 1 : length(F)-1
        if (force_to_right)
            % find displacements in direction of f, dx
            dx = a/2*(cos(theta2(i+1))-cos(theta2(i)));
            %find small amount of Wf
            Wf(i) = F(i+1)*dx;
            % Accumulate Wf
            Wf_sum = Wf_sum + Wf(i);
        end
    end
else
    Wf_sum = 0;
end

if (gravity)
    % Determine COM heihgts at initial SEP, state 1
    theta2_S = [1.62035622193624];
    thet33_S = [0.163640985321142];
    theta4_S = [1.87953152031875];
    
    y2_state1 = (a/2)*sin(theta2(1));
    y3_state1 = a*sin(theta2(1))+(b/2)*sin(theta3(1));
    y4_state1 = c/2*sin(theta4(1));
    
    % Determine COM heihgts at final SEP, state 2
    y2_state2 = (a/2)*sin(theta2(end));
    y3_state2 = a*sin(theta2(end))+(b/2)*sin(theta3(end));
    y4_state2 = c/2*sin(theta4(end));
    
    % Determine the GPE at initial and final SEPs
    Egrav_state2 = (y2_state2-y2_state1)*(abs(Fg2)) + (y3_state2-y3_state1)*abs(Fg3) + (y4_state2-y4_state1)*abs(Fg4);
    Egrav_state1 = 0;
else
    Egrav_state2 = 0;
    Egrav_state1 = 0;
end

% Find spring energy and work as a function of theta. First find net
% deflections at inital SEP states
phiB_state2 = theta3(end) - theta2(end) - thetaB;
phiC_state2 = theta4(end) - theta3(end) - thetaC;
phiB_state1 = theta3(1) - theta2(1) - thetaB;
phiC_state1 = theta4(1) - theta3(1) - thetaC;

% Calculate spring energy at both SEP states
Espring_state2 = .5*kB*phiB_state2^2 +.5*kC*phiC_state2^2;
Espring_state1 = .5*kB*phiB_state1^2 +.5*kC*phiC_state1^2;

% Calculate total energy at both states
E_state2 = Egrav_state2 + Espring_state2;
E_state1 = Egrav_state1 + Espring_state1;

% Print energy results
fprintf('\nE_state2 - E_state1 = %.4f\n',E_state2 - E_state1)
fprintf('Wf = %.4f\n',Wf_sum)
fprintf('E_state1 = %.4f\n',E_state1)
fprintf('E_state2 = %.4f\n',E_state2)
fprintf('Egrav_state1 = %.4f\n',Egrav_state1)
fprintf('Egrav_state2 = %.4f\n',Egrav_state2)
fprintf('Espring_state1 = %.4f\n',Espring_state1)
fprintf('Espring_state2 = %.4f\n',Espring_state2)

%% Determine the final configuration
% determine final angles
[theta3,theta4,g,h] = fourbar(theta2(end),a,b,c,d);

% determine location of pins at final configuration (m)
x(1) = 0;
y(1) = 0;
x(2) = a*cos(theta2(end));
y(2) = a*sin(theta2(end));
x(3) = x(2) + b*cos(theta3);
y(3) = y(2) + b*sin(theta3);
x(4) = d;
y(4) = 0;

% print results in command window
fprintf('\nSolution Angles:\n')
fprintf('theta2: %1.4f rad\ntheta3: %1.4f rad\ntheta4: %1.4f rad\n',theta2(end),theta3,theta4)

%% Plot results 
%test case 1 - sicilano
% x_sicilano = [0,0.766900797916280,-3.11930591468120,-2.00000000000000]+d;
% y_sicilano = [0,3.92579459041802,4.87310519785482,-8.88178419700125e-16];
% sicilano = plot(x_sicilano,y_sicilano,'r--','LineWidth',2)
% axis([-1.5 6.5 0 5])
% legend([constans,sicilano],'Constans approach','Sicilano approach','Location','SOUTHEAST')
% theta_text = sprintf('\\theta_{4} = \\theta_{1} = %.4f rad', theta4);
% text(2.75,2, theta_text, 'FontSize', 10);         % put equation on plo


%test case 2
x_sicilano_S = [0,-1.21541568510784,-6.14861902169833,-6]+d;
y_sicilano_S = [0,3.81087453380374,2.99631646966562,1.33226762955019e-15];
x_sicilano_S_G = [0,-1.75706862770346,-6.71077173724050,-6.00000000000000]+d;
y_sicilano_S_G = [0,3.59342591930602,2.91458462521508,-4.44089209850063e-16];

x_costans_S = [0,-0.148619061270426,4.78458427699186,6];
y_constans_S = [0,2.99631646770282,3.81087452171606,0];
x_costans_S_G = [0,-0.710771492787956,4.24293160905407,6]
y_constans_S_G = [0,2.91458468482904,3.59342603507275,0];

theta4_grav = [1.80999462932870];
theta4_no_grav = [1.62035629951925];

costans_S = plot(x_costans_S,y_constans_S,'k-','LineWidth',3);
costans_S_G = plot(x_costans_S_G,y_constans_S_G,'k--','LineWidth',3);
sicilano_S = plot(x_sicilano_S,y_sicilano_S,'m-','LineWidth',0.5);
sicilano_S_G = plot(x_sicilano_S_G,y_sicilano_S_G,'m--','LineWidth',0.5);



theta_grav_text = sprintf('\\theta_{4,Grav} = \\theta_{1,Grav} = %.4f rad', theta4_grav);
text(2,4.3, theta_grav_text, 'FontSize', 10);         % put equation on plo
theta_no_grav_text = sprintf('\\theta_{4,NoGrav} = \\theta_{1,NoGrav} = %.4f rad', theta4_no_grav);
text(2,4.75, theta_no_grav_text, 'FontSize', 10);         % put equation on plo

title_text = sprintf('Test Case with [a,b,c,d] = [%.f, %.f, %.f, %.f]',a,b,c,d);
xlabel('x-position (m)')
ylabel('y-position (m)')
axis equal
axis([-1.5 6.5 0 5])
legend([costans_S,costans_S_G,sicilano_S,sicilano_S_G],'Constans: No Grav.','Constans: Grav.','Sicilano: No Grav.','Sicilano: Grav.', 'Location','NORTHWEST')

title(title_text)
