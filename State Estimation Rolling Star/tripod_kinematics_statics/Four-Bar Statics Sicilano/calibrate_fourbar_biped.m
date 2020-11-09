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
% Now, span to our max seconds
max_t = 2;
times = 0:dt:max_t;

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













