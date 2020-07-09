function [simtimes, vertex_xy] = load_simRollingStar_dataset(fname)
%load_simRollingStar_dataset Load in CSV data from the DER simulations for
%all nodes of the robot.
%   
%   load_simRollingStar_dataset reads and parses CSV log files from the DER
%   simulations and organizes them into a dataset that's useful for
%   kinematics tests
%
%   Inputs:
%
%   fname = name of file to read. Assume for now it's in the same
%       directory as this function.
% 
%   Outputs:
%
%   simtimes = vector of the time points for each observation. Row
%       correlates with vert_xy.
%
%   vertex_xy = coordinates for each vertex, each timestep.
%
%   *TO-DO: might be easier to read in the "simulation options" from these
%   files, but Drew tried out some variations on textscan and readtable,
%   both seemed overkill and frustrating. We can hard-code the limb length,
%   etc.

% Hard-coded for our simulation: data starts at line 6, so counting from
% zero, 
simdata = csvread(fname, 6, 0);
% split into timestamps and vertex positions
simtimes = simdata(:,1);
vertex_xy = simdata(:,2:end);

end










