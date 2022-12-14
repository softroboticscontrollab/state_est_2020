% biped_statics.m
% Copyright 2020 Sam Alvares, Andrew Sabelhaus, and the Soft Machines Lab
% at Carnegie Mellon University

% Adapted from Sam's newton_four_bar_statics, this script tests a set of
% calibrated spring constants against the biped robot.

%% Prep the workspace
clc
clear all
close all

% We need some functions in other folders
addpath( genpath('../Tripod DER Simulation') );

%% Read in the DER data.

% Manually specify the filename. 
der_data_fname = '../Tripod DER Simulation/simBipedAllNodes_2020_10_23_163509.csv';
% Don't start at time zero, we get infinite radii / zero curvature since
% the limbs are flat
t_0 = 1.1;
% Now, span to our max seconds
max_t = 3.1;
[kappa, q, a, s, masses, times, vertices, feet] = get_biped_traj(der_data_fname, t_0, max_t);

%% Enter parameters for numeric solver
n = 4;      % number bars

% Pick out one timepoint for now.
timept_sec = 2 - t_0; % sec
% Assuming a 5 msec sampling time...
timept = round(200*timept_sec)+1;

% Curvatures:
kappa_t = cell2mat(kappa{timept})';
% bar lengths: it's "s" here except for...
a4 = a{timept}(4);
% masses specified manually
m1 = masses(1);
m2 = masses(2);
m3 = masses(3);
g=9.81;     % acceleration due to gravity (m/s/s)

% From our least-squares fit, hard-coded:
% Trying it with forced zero rest angle
k2 = 5.73995885893841;     % spring constant at joint 2 (m N/rad)
k3 = 5.76522370043459;     % spring constant at joint 3 (m N/rad)
q_2_bar = pi/2;
q_3_bar = pi/2;

% determine bar lengths
%K = [kappa1; kappa2; kappa3]; % create vector of curvature
K = kappa_t;
%L = [L1; L2; L3];             % create vector of arc lengths
L = s;
a = zeros(n-1,1);
for i = 1:n-1
    a(i) = barcalc(K(i),L(i)); % caluclate virtual bar length
end

%calculate initial spring lengths at t2i, the initial guess for t2
t1i=pi/2.3;  % initial guess of t (rad)
x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1i));    
alpha = asin((a4*sin(t1i))/x);
beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
% zeta = asin((x*sin(beta))/a(3));
zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
t2i = pi - alpha - beta;  % initial t2 (rad)
t3i = pi - zeta;          % initial t3 (rad)

% Sam's approach set the initial guess as the rest angle:
%t02 = pi -(t2i);          % unstretched spring length for t2 (rad)
%t03 = pi - (t3i);         % unstretched spring length for t3 (rad)
% Now we can use our fitted example:
t02 = q_2_bar;
t03 = q_3_bar;

%create a vector of parameters for the numeric solver
param = [s(1) s(2) s(3) a4 K(1) K(2) K(3) m1 m2 m3 g k2 k3 t02 t03]; %create a vector of all paramters

%% Run Newton Raphson numeric solver
R=@resid_four_bar_statics;      % residual function
dRdx=@dresid_four_bar_statics;  % slope of residual
tol=0.00001;                    % solution tolerance
maxIter=1000;                   % max iteratios to find solution
toggle=1;                       % 1 prints guesses

% Run numeric solver
[t1_sol,er_est,t1_guess]=func_newton(R,dRdx,t1i,tol,maxIter,toggle,param);

%% Forward Kinematics
% find the numeric value of t2 and t3 using the output of the numeric
% solver
x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1_sol));
alpha = asin((a4*sin(t1_sol))/x);
beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
t2_sol = pi - alpha - beta; % settling angle for joint 2 (rad)
t3_sol = pi - zeta;         % settling angle for joint 3 (rad)

% Determine forward kinamtics at the final pose
theta = [t1_sol;t2_sol; t3_sol];
zero = [0;0;0];
[T0_n,Tnm1_n] = fwdkinRISSnum(a, zero, zero, theta);

% Create corner vector at final position
x(1) = 0;
y(1) = 0;
x(2) = a(1)*cos(t1_sol);
y(2) = a(1)*sin(t1_sol);
x(3) = x(2) + a(2)*cos(t1_sol+t2_sol);
y(3) = y(2) + a(2)*sin(t1_sol+t2_sol);
x(4) = x(3) + a(3)*cos(t1_sol+t2_sol+t3_sol);
y(4) = y(3) + a(3)*sin(t1_sol+t2_sol+t3_sol);
x(5) = 0;
y(5) = 0;

% Determine forward kinamtics at the inital pose
theta_initial = [t1i;t2i;t3i];
zero = [0;0;0];
[T0_n_initial,Tnm1_n_initial] = fwdkinRISSnum(a, zero, zero, theta_initial);

% Create corner vector for initial pose
xin(1) = 0;
yin(1) = 0;
xin(2) = a(1)*cos(t1i);
yin(2) = a(1)*sin(t1i);
xin(3) = xin(2) + a(2)*cos(t1i+t2i);
yin(3) = yin(2) + a(2)*sin(t1i+t2i);
xin(4) = xin(3) + a(3)*cos(t1i+t2i+t3i);
yin(4) = yin(3) + a(3)*sin(t1i++t2i+t3i);
xin(5) = 0;
yin(5) = 0;

%% Find COM Locations
% Determine the location of the COM at the final pose
P0Gx = zeros(n-1,1);
P0Gy = zeros(n-1,1);
for i = 1:n-1
    phi = L(i)*K(i)*.5;
    A =(2*sin(phi))/(L(i)*K(i)^2);
    B = cos(phi)/K(i);
    C = A-B;
    if i==1
        P0Gx(1) = a(1)/2*T0_n{1}(1,1)+C*T0_n{1}(1,2);
        P0Gy(1) = a(1)/2*T0_n{1}(2,1)+C*T0_n{1}(2,2);
    else
        P0Gx(i) = T0_n{i-1}(1,4)+a(i)/2*T0_n{i}(1,1)+C*T0_n{i}(1,2);
        P0Gy(i) = T0_n{i-1}(2,4)+a(i)/2*T0_n{i}(2,1)+C*T0_n{i}(2,2);
    end
end

% Determine the location of the COM at the initial pose
P0Gx_initial = zeros(n-1,1);
P0Gy_initial = zeros(n-1,1);
for i = 1:n-1
    phi = L(i)*K(i)*.5;
    A =(2*sin(phi))/(L(i)*K(i)^2);
    B = cos(phi)/K(i);
    C = A-B;
    if i==1
        P0Gx_initial(1) = a(1)/2*T0_n_initial{1}(1,1)+C*T0_n_initial{1}(1,2);
        P0Gy_initial(1) = a(1)/2*T0_n_initial{1}(2,1)+C*T0_n_initial{1}(2,2);
    else
        P0Gx_initial(i) = T0_n_initial{i-1}(1,4)+a(i)/2*T0_n_initial{i}(1,1)+C*T0_n_initial{i}(1,2);
        P0Gy_initial(i) = T0_n_initial{i-1}(2,4)+a(i)/2*T0_n_initial{i}(2,1)+C*T0_n_initial{i}(2,2);
    end
end

%% Plot results
% Plot initial pose
hold on 
initial_corners = plot(xin,yin,'--k');
axis equal
%axis([-7 2 -2 7])

% Plot COM at inital pose 
hold on
initial_COM = plot(P0Gx_initial,P0Gy_initial,'kx');

% Plot final pose
hold on
final_corners = plot(x,y,'r');

% Plot COM at final pose 
hold on
final_COM = plot(P0Gx,P0Gy,'rx');

% %plot limb verticies from data structure
% hold on
% plot(vertices{timept}(1,:), vertices{timept}(2,:),'gx')

% % plot the feet verticies from data structure
% for i = 1:num_feet
%     hold on
%     plot(feet{timept}{i}(1,:), feet{timept}{i}(2,:),'gx')
% end
% axis([-.02 0.02 0 .03])
% axis equal

% Create legend
title('Biped Statics Example');
legend('Initial Guess','Guess COM','Soln. Pose','Soln. COM')
xlabel('x-position (m)')
ylabel('y-position (m)')


% Determine unactuated pose 
% xi(1) = 0;
% yi(1) = 0;
% xi(2) = L(1)*cos(pi/2);
% yi(2) = L(1)*sin(pi/2);
% xi(3) = xi(2) + L(2)*cos(pi);
% yi(3) = yi(2) + L(2)*sin(pi);
% xi(4) = -a4;
% yi(4) = 0;
% xi(5) = 0;
% yi(5) = 0;

% determine inital COM pos
% P0Gxi = [0;-L2/2;-L2];
% P0Gyi = [L1/2;L1;L1/2];

% Plot non-actuated pose
% hold on
% plot(xi,yi,'--k')

%Plot COM position at unactuated position
% hold on
% plot(P0Gxi,P0Gyi,'kx')

%% Determine joint torques 
% % Gravity
% go = [ g*m1*(sin(t1_sol)*(cos((L1*kappa1)/2)/kappa1 - (2*sin((L1*kappa1)/2))/(L1*kappa1^2)) + (sin((L1*kappa1)/2)*cos(t1_sol))/kappa1) - g*m3*((cos(t1_sol + t2_sol + t3_sol)*sin((L3*kappa3)/2))/kappa3 - (2*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1 - (2*sin((L3*kappa3)/2)*cos(t1_sol + t2_sol)*cos(t3_sol))/kappa3 + (2*sin((L3*kappa3)/2)*sin(t1_sol + t2_sol)*sin(t3_sol))/kappa3 - (2*sin((L2*kappa2)/2)*cos(t1_sol)*cos(t2_sol))/kappa2 + (2*sin((L2*kappa2)/2)*sin(t1_sol)*sin(t2_sol))/kappa2 + (sin(t1_sol + t2_sol + t3_sol)*(2*sin((L3*kappa3)/2) - L3*kappa3*cos((L3*kappa3)/2)))/(L3*kappa3^2)) + g*m2*((sin((L2*kappa2)/2)*cos(t1_sol + t2_sol))/kappa2 + (2*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1 - (sin(t1_sol + t2_sol)*(2*sin((L2*kappa2)/2) - L2*kappa2*cos((L2*kappa2)/2)))/(L2*kappa2^2));
%                                                                                                                                                                                                                                                                                                           g*m2*(sin(t1_sol + t2_sol)*(cos((L2*kappa2)/2)/kappa2 - (2*sin((L2*kappa2)/2))/(L2*kappa2^2)) + (sin((L2*kappa2)/2)*cos(t1_sol + t2_sol))/kappa2) + (g*m3*(kappa2*cos(t1_sol + t2_sol + t3_sol + (L3*kappa3)/2) - kappa2*cos(t1_sol + t2_sol + t3_sol - (L3*kappa3)/2) - L3*kappa3^2*sin(t1_sol + t2_sol - (L2*kappa2)/2) + L3*kappa3^2*sin(t1_sol + t2_sol + (L2*kappa2)/2) + L3*kappa2*kappa3*sin(t1_sol + t2_sol + t3_sol + (L3*kappa3)/2)))/(L3*kappa2*kappa3^2);
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 (g*m3*(cos(t1_sol + t2_sol + t3_sol + (L3*kappa3)/2) - cos(t1_sol + t2_sol + t3_sol - (L3*kappa3)/2) + L3*kappa3*sin(t1_sol + t2_sol + t3_sol + (L3*kappa3)/2)))/(L3*kappa3^2)];
% ga = g*m1*(sin(t1_sol)*(cos((L1*kappa1)/2)/kappa1 - (2*sin((L1*kappa1)/2))/(L1*kappa1^2)) + (sin((L1*kappa1)/2)*cos(t1_sol))/kappa1) - g*m2*((sin((L2*kappa2)/2)*cos(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t1_sol + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)))))/kappa2 - (2*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1 + (sin(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t1_sol + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*(2*sin((L2*kappa2)/2) - L2*kappa2*cos((L2*kappa2)/2)))/(L2*kappa2^2)) - ((a4*(4*cos(t1_sol) - 4*cos((L1*kappa1)/2)^2*cos(t1_sol) + 2*a4*kappa1*sin((L1*kappa1)/2) + a4^2*kappa1^2*cos(t1_sol) + 2*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol)^2))/(kappa1^2*(1 - (a4^2*kappa1^2*sin(t1_sol)^2)/(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2))^(1/2)*((a4^2*kappa1^2 - 4*cos((L1*kappa1)/2)^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol) + 4)/kappa1^2)^(3/2)) + (a4*sin((L1*kappa1)/2)*sin(t1_sol)*(4*kappa1^2*kappa2^2 - 4*kappa1^2*kappa3^2 + 4*kappa2^2*kappa3^2 - 4*kappa2^2*kappa3^2*cos((L1*kappa1)/2)^2 + 4*kappa1^2*kappa3^2*cos((L2*kappa2)/2)^2 - 4*kappa1^2*kappa2^2*cos((L3*kappa3)/2)^2 + a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2)*cos(t1_sol)))/(2*kappa1^3*kappa2*kappa3^2*sin((L2*kappa2)/2)*((a4^2*kappa1^2 - 4*cos((L1*kappa1)/2)^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol) + 4)/kappa1^2)^(3/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2)))*(g*m2*(sin(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t1_sol + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*(cos((L2*kappa2)/2)/kappa2 - (2*sin((L2*kappa2)/2))/(L2*kappa2^2)) - (sin((L2*kappa2)/2)*cos(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t1_sol + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)))))/kappa2) - (g*m3*(kappa2*cos(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) + (L3*kappa3)/2) - kappa2*cos(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) - (L3*kappa3)/2) + L3*kappa3^2*sin(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t1_sol + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) + (L2*kappa2)/2) + L3*kappa3^2*sin(t1_sol - asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) + (L2*kappa2)/2) + L3*kappa2*kappa3*sin(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) - (L3*kappa3)/2)))/(L3*kappa2*kappa3^2)) - g*m3*((cos(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*sin((L3*kappa3)/2))/kappa3 - (2*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1 + sin(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t1_sol + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2) - (sin(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*(2*sin((L3*kappa3)/2) - L3*kappa3*cos((L3*kappa3)/2)))/(L3*kappa3^2) + (2*cos(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*sin((L2*kappa2)/2)*cos(t1_sol))/kappa2 + (2*sin(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*sin((L2*kappa2)/2)*sin(t1_sol))/kappa2 - (2*sin((L3*kappa3)/2)*cos(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t1_sol + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))))*((kappa3^2*((kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)) - 1)*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2))/(4*kappa1^2*sin((L3*kappa3)/2)^2) + 1)^(1/2))/kappa3) + (a4*g*m3*sin((L1*kappa1)/2)*sin(t1_sol)*(cos(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) + (L3*kappa3)/2) - cos(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) - (L3*kappa3)/2) + L3*kappa3*sin(asin((kappa3*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2))/(2*sin((L3*kappa3)/2))) - t1_sol + asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))) - (L3*kappa3)/2))*(a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*cos(t1_sol)*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2) - 4*kappa1^2*kappa2^2*sin((L3*kappa3)/2)^2 - 4*kappa1^2*kappa3^2*sin((L2*kappa2)/2)^2 + 4*kappa2^2*kappa3^2*sin((L1*kappa1)/2)^2))/(L3*kappa1^3*kappa3^3*sin((L2*kappa2)/2)^2*sin((L3*kappa3)/2)*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2)*((a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*cos(t1_sol)*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2) - 4*kappa1^2*kappa2^2*sin((L3*kappa3)/2)^2 - 4*kappa1^2*kappa3^2*sin((L2*kappa2)/2)^2 + 4*kappa2^2*kappa3^2*sin((L1*kappa1)/2)^2)^2/(kappa1^4*kappa2^2*kappa3^2*sin((L2*kappa2)/2)^2*sin((L3*kappa3)/2)^2))^(1/2));
% % Spring forces 
% to = [                   0
%          -k2*(t02 + t2_sol - pi)
%          -k3*(t03 + t3_sol - pi)];
% ta = (a4*k3*sin((L1*kappa1)/2)*sin(t1_sol)*(t03 - asin((kappa3*((4*sin((L1*kappa1)/2)^2 + a4^2*kappa1^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(4*sin((L1*kappa1)/2)^2 + a4^2*kappa1^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol))))^(1/2))/(2*sin((L3*kappa3)/2))))*(4*kappa2^2*kappa3^2*sin((L1*kappa1)/2)^2 - 4*kappa1^2*kappa3^2*sin((L2*kappa2)/2)^2 - 4*kappa1^2*kappa2^2*sin((L3*kappa3)/2)^2 + a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2)*cos(t1_sol)))/(kappa1^3*kappa3*sin((L2*kappa2)/2)^2*sin((L3*kappa3)/2)*((a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2)*((a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*cos(t1_sol)*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2) - 4*kappa1^2*kappa2^2*sin((L3*kappa3)/2)^2 - 4*kappa1^2*kappa3^2*sin((L2*kappa2)/2)^2 + 4*kappa2^2*kappa3^2*sin((L1*kappa1)/2)^2)^2/(kappa1^4*kappa2^2*kappa3^2*sin((L2*kappa2)/2)^2*sin((L3*kappa3)/2)^2))^(1/2)) - k2*((a4*(4*cos(t1_sol) - 4*cos((L1*kappa1)/2)^2*cos(t1_sol) + 2*a4*kappa1*sin((L1*kappa1)/2) + a4^2*kappa1^2*cos(t1_sol) + 2*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol)^2))/(kappa1^2*(1 - (a4^2*kappa1^2*sin(t1_sol)^2)/(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2))^(1/2)*((a4^2*kappa1^2 - 4*cos((L1*kappa1)/2)^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol) + 4)/kappa1^2)^(3/2)) + (a4*sin((L1*kappa1)/2)*sin(t1_sol)*(4*kappa1^2*kappa2^2 - 4*kappa1^2*kappa3^2 + 4*kappa2^2*kappa3^2 - 4*kappa2^2*kappa3^2*cos((L1*kappa1)/2)^2 + 4*kappa1^2*kappa3^2*cos((L2*kappa2)/2)^2 - 4*kappa1^2*kappa2^2*cos((L3*kappa3)/2)^2 + a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2)*cos(t1_sol)))/(2*kappa1^3*kappa2*kappa3^2*sin((L2*kappa2)/2)*((a4^2*kappa1^2 - 4*cos((L1*kappa1)/2)^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1_sol) + 4)/kappa1^2)^(3/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1_sol)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2)))*(asin((a4*sin(t1_sol))/((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2)) - t02 + acos((kappa2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1))/(4*sin((L2*kappa2)/2)*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1_sol))/kappa1)^(1/2))));
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
%% Verify solution with energy
% Determine initial PE with unstretched spring lengths
% t1_guess = t1_guess(1:length(t1_guess)-1);
% for i = 1:length(t1_guess)
%     %find the numeric value of t2 and t3
%     x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1_guess(i)));
%     alpha = asin((a4*sin(t1_guess(i)))/x);
%     beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
%     % zeta = asin((x*sin(beta))/a(3));
%     zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
%     t2_guess = pi - alpha - beta;
%     t3_guess = pi - zeta;
%     % t3_sol = atan2(s_zeta,c_zeta);
%     
%     %determine forward kinamtics
%     theta = [t1_guess(i);t2_guess; t3_guess];
%     zero = [0;0;0];
%     [T0_n,Tnm1_n] = fwdkinRISSnum(a, zero, zero, theta);
%     
%     for j = 1:n-1
%         phi = L(j)*K(j)*.5;
%         A =(2*sin(phi))/(L(j)*K(j)^2);
%         B = cos(phi)/K(j);
%         C = A-B;
%         if j==1
%             P0Gx(1) = a(1)/2*T0_n{1}(1,1)+C*T0_n{1}(1,2);
%             P0Gy(1) = a(1)/2*T0_n{1}(2,1)+C*T0_n{1}(2,2);
%         else
%             P0Gx(j) = T0_n{j-1}(1,4)+a(j)/2*T0_n{j}(1,1)+C*T0_n{j}(1,2);
%             P0Gy(j) = T0_n{j-1}(2,4)+a(j)/2*T0_n{j}(2,1)+C*T0_n{j}(2,2);
%         end
%     end
%     
%     Ugrav(i) = g*(m1*P0Gy(1)+m2*P0Gy(2)+m3*P0Gy(3));
%     Uspring(i) = .5*k2*(pi-t2_guess-t02)^2+.5*(pi-t3_guess-t03)^2;
%     Utot(i) = Ugrav(i) + Uspring(i);
% end
% Ugrav
% Uspring
% Utot
% 
% x_axis = 1:length(t1_guess);
% figure
% subplot(3,1,1);
% plot(x_axis,Ugrav,'b-')
% ylabel('Ugrav')
% hold on
% subplot(3,1,2)
% plot(x_axis,Uspring,'k-')
% ylabel('Uspring')
% hold on
% subplot(3,1,3)
% plot(x_axis,Utot,'r-')
% ylabel('Utot')
% 
% t2_constans = 2*pi-(pi-t1_sol)-(pi-t2_sol)-(pi-t3_sol);




