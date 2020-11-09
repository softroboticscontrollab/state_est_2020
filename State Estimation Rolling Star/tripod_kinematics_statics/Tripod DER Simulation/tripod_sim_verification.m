%% tripod_sim_verification.m
% Sam Alvares and Drew Sabelhaus
% 11/2/2020
% This script will read simulation data (from the discrete elastic rods
% simulation) of the rolling star robot and create a plot to verify that
% tripod_sim_parse.m is functioning properly

%% Prep the workspace
clear all;
close all;
clc;

%% Define constants
num_limbs = 3;  % number of limbs
num_feet = 2;   % number of feet
   
%% Read in data and parse
% Manually specify the filename. We'll assume the file is in the same
% directory as the m-file that loads it.
der_data_fname = 'simBipedAllNodes_2020_10_23_163509.csv';
simdata = csvread(der_data_fname, 5, 0);
simtime = simdata(:,1); % sec
    
%%  Parse data from simulation at a particular time
time = 2;       % time of simulation simulation to analyze
[vertex_xy,feet_data,limbs_data] = tripod_sim_parse(simdata,time);

%% Plot the outputs 
% Plot all verticies
figure(1)
plot(vertex_xy(1,:),vertex_xy(2,:),'k.')
axis([-.02 0.02 0 .03])
axis equal

% Plot limb verticies and foot verticies
figure(2)
% plot limb verticies from data structure
for i = 1:3
    hold on
    plot(limbs_data{i}(1,:),limbs_data{i}(2,:),'rx')
end
% plot the feet verticies from data structure
for i = 1:num_feet
    hold on
    plot(feet_data{i}(1,:),feet_data{i}(2,:),'kx')
end
axis([-.02 0.02 0 .03])
axis equal
