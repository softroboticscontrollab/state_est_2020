% star_kin_sim_verification.m
% (c) 2020 Andrew P. Sabelhaus, Sam Alvares, and the Soft Machines Lab at
% Carnegie Mellon University

% This script will read simulation data (from the discrete elastic rods
% simulation) of the rolling star robot, and verify the reduced-order
% kinematics model against the simulation.

%% Prep the workspace

clear all;
close all;
clc;

% The function to load from the simulator's CSV file is in a subdirectory
addpath( genpath('simRollingStarData') );

%% Load data from the DER simulator

% Manually specify the filename. We'll assume the file is in the same
% directory as the m-file that loads it.

der_data_fname = 'simRollingStarAllNodes_2020_07_09_143539.csv';
[simtimes, vertex_xy] = load_simRollingStar_dataset(der_data_fname);

% The vertex coordinates are stored as [x1, y1, x2, y2, .... xN, yN], but
% we'll need to separate out the x and y values. Take every other column:
vertex_x = vertex_xy(:, 1:2:end);
vertex_y = vertex_xy(:, 2:2:end);

% Hard-code some of the simulation parameters from the text file.
num_limbs = 7;
num_v_per_circ = 13; % number of vertices per "circular part" of a limb
num_v_per_flat = 5; % number of vertices per "flat part" of a limb.
% ...note that total number of vertices, i.e. columns in vertex_x or vertex_y, is:
% num_limbs * (num_v_per_circ + num_v_per_flat).


%% Do something as an example

% Plot all nodes of the robot at the first timestep.
figure;
hold on;
plot(vertex_x(1,:), vertex_y(1,:), 'k.');
title('All nodes at time zero');

% since the data is a bit confusing in terms of "circular" versus "flat"
% sections, here's an example of each.
timezero_circ_x = vertex_x(1, 1:(num_limbs * num_v_per_circ));
timezero_circ_y = vertex_y(1, 1:(num_limbs * num_v_per_circ));
timezero_flat_x = vertex_x(1, (num_limbs * num_v_per_circ):end);
timezero_flat_y = vertex_y(1, (num_limbs * num_v_per_circ):end);

figure;
hold on;
plot(timezero_circ_x, timezero_circ_y, 'r.');
title('Circular sections of the robot only');

figure;
hold on;
plot(timezero_flat_x, timezero_flat_y, 'b.');
title('Flat sections (sticking out past the circular sections)');






