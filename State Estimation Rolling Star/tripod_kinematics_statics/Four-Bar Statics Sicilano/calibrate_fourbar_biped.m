% calibrate_fourbar_biped.m
% Copyright 2020 Andrew P. Sabelhaus and the Soft Machines Lab at Carnegie
% Mellon University

% This script calculates calibrated spring constants k_i and angular
% spring rest angles \bar q_i for the four-bar approximation to the
% biped robot from DER simulations.

% Approach outline:
% 1) Read vertex data from simulation run. Calculate curvatures kappa_i and
% virtual joint angles q_i for all timepoints t=1...T. Also extract arc
% length.
% 2) Formulate a linear fit problem. TO-DO: HOW TO GET THIS TO BE
% ONLY TWO SPRINGS?
% 3) Solve via least squares / pseudoinverse. Calculate \bar q_i from the
% rest joint torques
% 4) Save the results

%% Prep the workspace

clear all;
close all;
clc;

% We need some functions in other folders
addpath( genpath('../Tripod DER Simulation') );

%% 1) Read vertex data and calculate reduced state

% Manually specify the filename. 
der_data_fname = '../Tripod DER Simulation/simBipedAllNodes_2020_10_23_163509.csv';
simdata = csvread(der_data_fname, 5, 0);
% We want all the timepoints. The function takes in time in sec, so we need
% to specify a dt to span out everything we want
dt = 0.005;
% Don't start at time zero, we get infinite radii / zero curvature since
% the limbs are flat
t_0 = 1.1;
% Now, span to our max seconds
max_t = 3.1;
times = t_0:dt:max_t;

% And extract vertices, etc. from each datapoint
vertices = {};
feet = {};
limbs = {};
for t = 1:size(times, 2)
    % This timestep
    [vertex_xy_t, feet_data_t, limbs_data_t] = tripod_sim_parse(simdata,times(t));
    % Add to our running lists
    vertices{end+1} = vertex_xy_t;
    feet{end+1} = feet_data_t;
    limbs{end+1} = limbs_data_t;
end

% We want to flip around the limb and foot numbering here. The rightmost leg should
% be limb 1.
% Here's an ugly manual way.
for t = 1:size(limbs,2)
    limb_1_from_3 = limbs{t}{3};
    limbs{t}{3} = limbs{t}{1};
    limbs{t}{1} = limb_1_from_3;
    % and the feet
    foot_1_from_2 = feet{t}{2};
    feet{t}{2} = feet{t}{1};
    feet{t}{1} = foot_1_from_2;
end

% Calculate our required variables. First, some constants: 
n = 4;          % number bars
s_i = 0.022;    % arc length, all limbs (m)
s = repmat(s_i, n-1, 1); % patterned out per limb
g = 9.81;       % acceleration due to gravity (m/s/s)
n_per_limb = 13;% number of nodes for one limb

% Mass of one limb: given a density, rod radius, and arc length,
density = 1912.00;
radius = 2e-3;
% the mass per limb is
m_i = pi * radius*2 * s_i * density; % in kg.
% patterned out
masses = repmat(m_i, n-1, 1);

% Curvatures: will have the same format as limbs{},
% which is limbs{t}{limb_num} is a 2 x n_per_limb matrix.
kappa = {};

% calculate circular fit for each limb, each timestep. We have n-1 = 3 limbs.
for t = 1:size(times, 2)
    % each of our three limbs...
    for j = 1:n-1
        [xc,yc,Re] = circfit(limbs{t}{j}(1,:), limbs{t}{j}(2,:)); % create a curve fit for each of the limbs
        kappa{t}{j} = 1/Re; % record curvatures
        %th = linspace(0,2*pi,100)'; % create a vector of angles to plot circle
        % calcuale circle from curve fit
        %xe = Re*cos(th)+xc;
        %ye = Re*sin(th)+yc;
        %%plot curve fit circles
        %figure(2)
        %hold on
        %axis equal
        %plot([xe;xe(1)],[ye;ye(1)])
    end
end

% Get the virtual angles between bars and the virtual bar lengths
% Per Sam's code, a(1) is right-side bar, a(2) is top, a(3) is left, a(4)
% is the ground.
% And, q_i go bottom right, top right, top left, bottom left.
q = {};
a = {};

for t = 1:size(times, 2)
    % The corners here are
    corners_t = biped_corners(limbs{t}, feet{t});
    % Verify:
    %plot(corners_t(1,:), corners_t(2,:));
    % q_t = [q1_t; q2_t, q3_t, q4_t]
    q{t} = zeros(n, 1);
    q{t}(1) = pi - int_ang(corners_t(:,end), corners_t(:,1), corners_t(:,2));
    % the rest
    for m = 2:n-1
        q{t}(m) = pi - int_ang(corners_t(:,m-1), corners_t(:,m), corners_t(:,m+1));
    end
    q{t}(n) = pi - int_ang(corners_t(:,end-1), corners_t(:,end), corners_t(:,1));
    % Then the virtual bar lengths
    % For all limbs at this timestep
    kappa_t = cell2mat(kappa{t});
    % the virtual bar lengths
    a_t = zeros(n-1,1);
    % per limb bar length...
    for j = 1:n-1
        a_t(j) = barcalc(kappa_t(j), s(j));
    end
    % finally, a4 is the distance between the two corners.
    a_t4 = norm(corners_t(:,1) - corners_t(:,4));
    a_t(4) = a_t4;
    % save
    a{end+1} = a_t;
end

%% 2) Formulate the linear regression problem.

% Right-hand side:
% ga_t is a scalar,
% Ga is a T x 1 column vector of ga_t, t=1...T timesteps
T = size(times,2);
Ga = zeros(T,1);

% Left-hand side: 

% B_t needs to be a (n-1) x 2m matrix, with the first row zeros because
% Sam's code doesn't use a spring at q1 and also only includes three angles
% {q1, q2, q3}, no need for q4.
% The (n-1) is number of joint angles (unconstrained), and the m is number 
% of springs because we have both a k_i and tau_bar_i per spring.
num_springs = 2;
B_t = zeros(n-1, 2*num_springs);

% When applying Y, we'll get B_tilde_t, which is 1 x 2m, then H, which is
% diag(B_tilde_1, ..., B_tilde_T), size T*(2mT)
%H = zeros(T, 2*num_springs*T);
% ...nevermind, we can use MATLAB's functions here instead of manually
% indexing.
H = [];

for t=1:T
    % Build up the parameters list required for the symbolically-solved
    % functions.
    % L1=param(1);
    % L2=param(2);
    % L3=param(3);
    % a4=param(4);
    % kappa1 = param(5);
    % kappa2 = param(6);
    % kappa3 = param(7);
    % m1=param(8);
    % m2=param(9);
    % m3=param(10);
    % g=param(11);
    % t1=param(12);
    
    % The kappas are stored weird but whatever. It's a 1 x T cell array
    param_t = [ s; 
                a{t}(4);
                cell2mat(kappa{1,t})';
                masses;
                g;
                q{t}(1)];
    % Right-hand side
    Ga(t) = ga_fourbar(param_t);
    % Left-hand side
    Y_t = closure_Y_fourbar(param_t);
    B_t = zeros(n-1, 2*num_springs);
    % Insert the angles, bottom-left corner
    B_t(2:end, 1:2) = diag([q{t}(2); q{t}(2)]);
    % The rest torques are bottom-right corner
    B_t(2:end, 3:4) = -eye(2);
    % disp('B_t is');
    % B_t
    % Apply the loop closure constraint
    B_tilde_t = Y_t' * B_t;
    % ...which is 1 x 4, and can be stacked into...
    H = blkdiag(H, B_tilde_t);
end

%% 3) Solve via least squares.

% This should be as simple as...
x = pinv(H) * Ga;
% x = [k2, k3, tau_bar2, tau_bar3]
% Pick out the rest angles from the rest torques as q_bar_i = tau_i / k_i
q_bar_2 = x(3) / x(1);
q_bar_3 = x(4) / x(2);

q_bar = [q_bar_2; q_bar_3];
% In degrees, for reference:
disp('Fitted spring constants are:');
x(1:2)
disp('Fitted rest angles for springs are:');
q_bar * 180/pi












