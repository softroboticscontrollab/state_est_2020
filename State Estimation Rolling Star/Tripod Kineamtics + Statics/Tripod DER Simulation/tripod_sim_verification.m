%% tripod_sim_verification.m
% This script will read simulation data (from the discrete elastic rods
% simulation) of the rolling star robot.

%% Prep the workspace
clear all;
close all;
clc;

%% Define constants
time = 0; % simulation time
num_limbs = 3; %number of limbs
num_v_per_lim = 13; % number of vertices per "circular part" of a limb
num_v_per_foot = 2; % number of vertices per "flat part" of a limb.
% ...note that total number of vertices, i.e. columns in vertex_x or vertex_y, is:
% num_limbs * (num_v_per_circ + num_v_per_flat).
   
%% Load data from the DER simulator
% The function to load from the simulator's CSV file is in a subdirectory
addpath(genpath('simRollingStarData') );

% Manually specify the filename. We'll assume the file is in the same
% directory as the m-file that loads it.
der_data_fname = 'simBipedAllNodes_2020_10_23_163509.csv';
simdata = csvread(der_data_fname, 5, 0);

simtime = simdata(:,1);
% 
% for i = 1:length(simtime)
%     sim_xy{i} = zeros(2,126);
%     sim_COM_xy{i} = zeros(2,1);
%     model_COM_xy{i} = zeros(2,1);
%     curve_xy{i} = zeros(2,7);
%     straight_xy{i} = zeros(2,7);
% end
% 
% for q = 1:length(simtime)
%     time = simtime(q);
    
%%  Parse data from simulation at a particular time
[vertex_xy,curve_data,straight_data,circ_tips,tips] = rollingStar_time(simdata,time);