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
% Don't start at time zero, we get infinite radii / zero curvature since
% the limbs are flat
t_0 = 1.995;
% Now, span to our max seconds
max_t = 2.005;
[kappa, q, a, s, masses, times, ~, ~] = get_biped_traj(der_data_fname, t_0, max_t);

% There are four virtual bars in the biped
n=4;
g = 9.81;       % acceleration due to gravity (m/s/s)

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
%B_t = zeros(n-1, 2*num_springs);
% Trying one with just the spring constants.
B_t = zeros(n-1, num_springs);

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
    % Trying without the torques
    %B_t = zeros(n-1, 2*num_springs);
    B_t = zeros(n-1, num_springs);
    % Insert the angles, bottom-left corner
    B_t(2:end, 1:2) = diag([q{t}(2); q{t}(2)]);
    % The rest torques are bottom-right corner
    % ...trying without the torques.
    %B_t(2:end, 3:4) = -eye(2);
    % disp('B_t is');
    % B_t
    % Apply the loop closure constraint
    B_tilde_t = Y_t' * B_t;
    % ...which is 1 x 4, and can be stacked into...
    H = blkdiag(H, B_tilde_t);
end

% Finally, last bit: we need to add a little something so we're solving for
% one set of constants, not T-many of them.

%I_bar = repmat(eye(num_springs*2), T, 1);
I_bar = repmat(eye(num_springs), T, 1);
% Then our lefthandside is Wx, where
W = H * I_bar;

%% 3) Solve via least squares.

% This should be as simple as...
x = pinv(W) * Ga;
% x = [k2, k3, tau_bar2, tau_bar3]
% Pick out the rest angles from the rest torques as q_bar_i = tau_i / k_i
%q_bar_2 = x(3) / x(1);
%q_bar_3 = x(4) / x(2);

%q_bar = [q_bar_2; q_bar_3];
% In degrees, for reference:
disp('Fitted spring constants are:');
x(1:2)
% disp('Fitted rest angles for springs in radians:');
% q_bar
% disp('Fitted rest angles for springs are, in degrees:');
% q_bar * 180/pi












