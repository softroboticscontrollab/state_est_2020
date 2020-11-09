function [vertex_xy,feet_data,limbs_data] = tripod_sim_parse(simdata,time)
% Sam Alvares and Drew Sabelhaus
% 11/2/2020
% tripod_sim_parse.m inputs the data from a tripod DER simulation and a
% simulation time. The function parses the simulation data into useful
% vectors and data structures
%
%   Inputs:
%
%   simdata = simulation data from cvs file.
%
%   time = point in simulation where we want to inspect the robot
% 
%   Outputs:
%
%   vertex_xy = 2 row matrix with all the x and y coordinates of all the
%   verticies (x in row 1, y in row 2)
%
%   feet_data = struct with 2x2 matricies of x and y coordinates of the feet 
%
%   limbs_data = struct with 13 or 14x2 matricies with 2x2 matricies of x 
%   and y coordinates of the limbs
%

% calulcate the row of the spreadsheet that is closest to the time input
row = round(200*time)+1;
simdata = simdata(row,:);

% Hard-code some of the simulation parameters from the text file.
num_limbs = 3;  % number of limbs
num_feet = 2;   % number of feet
num_v_per_limb = 13; % number of vertices per "circular part" of a limb
num_v_per_foot = 2; % number of vertices per "flat part" of a limb.
% ...note that total number of vertices, i.e. columns in vertex_x or vertex_y, is:
% num_limbs * (num_v_per_limb ) + num_feed * (num_v_per_foot).

%create a matrix excluding time stamp data
vertex_all_data = simdata(2:end);

% The vertex coordinates are stored as [x1, y1, x2, y2, .... xN, yN], but
% we'll need to separate out the x and y values. Take every other column:
vertex_x = vertex_all_data(1:2:end);
vertex_y = vertex_all_data(2:2:end);

% create a matrix of the x and y coordinates of all verticies. Row 1 is the
% x coordinate and row 2 is the y coordinate
vertex_xy = [vertex_x;vertex_y];

% Create vectors of the x and y positions of the limbs excluding the feet
limbs_x = vertex_x(1:(num_limbs * num_v_per_limb+1)); %x-coordinates for the limbs
limbs_y = vertex_y(1:(num_limbs * num_v_per_limb+1)); %y-coordinates for the limbs

% Create vectors of the x and y positions of the feet excluding the limbs
feet_x = vertex_x((num_limbs * num_v_per_limb)+2:end); %x-coordinates for the feet
feet_y = vertex_y((num_limbs * num_v_per_limb)+2:end); %y-coordinates for the feet

%create data structure for limbs
for i = 1:num_limbs    
    if i == 3 %the third limb has 14 data points
        limb_x_vec = limbs_x((i-1)*num_v_per_limb+1:(i)*num_v_per_limb+1);
        limb_y_vec = limbs_y((i-1)*num_v_per_limb+1:(i)*num_v_per_limb+1);
    else %all other limbs have 13 data points
        limb_x_vec = limbs_x((i-1)*num_v_per_limb+1:(i)*num_v_per_limb);
        limb_y_vec = limbs_y((i-1)*num_v_per_limb+1:(i)*num_v_per_limb);
    end    
    %create data structure of x and y coordinates for each limb
    limbs_data{i} = [limb_x_vec; limb_y_vec];
end

%create data structure for the feet
for i = 1:num_feet
    foot_x_vec = feet_x((i-1)*num_v_per_foot+1:(i)*num_v_per_foot); %x-coordinates for the curved sections
    foot_y_vec = feet_y((i-1)*num_v_per_foot+1:(i)*num_v_per_foot); %y-coordinates for the curved sections
    feet_data{i} = [foot_x_vec; foot_y_vec];
end

