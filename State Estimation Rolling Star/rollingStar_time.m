function [vertex_xy,curve_data,stright_data,curve_tips,tips] = rollingStar_time(simdata,time)
%,stright_data,circ_tips,tips
%rollingStar_row inputs the output data from load_simRollingStar_dataset
%and a specific time.  The function outputs a data structure with the x and
%y pos of the verticies separated by arm
%
%   Inputs:
%
%   simdata = parsed data from cvs file.
%   time = point in simulation where we want to inspect the robot
% 
%   Outputs:
%
%   curve_data = struct with 14x2 matricies with the xy verteixes of the
%   curved section of the arms. Both the first and 14th data points are
%   tips
%
%   stright_data = struct with 14x2 matricies with the xy verteixes of the
%   stright section of the arms
%

% calulcate the row of the spreadsheet that is closest to the time input
row = round(200*time+6);
simdata = simdata(row,:);

% Hard-code some of the simulation parameters from the text file.
num_limbs = 7;
num_v_per_circ = 13; % number of vertices per "circular part" of a limb
num_v_per_flat = 5; % number of vertices per "flat part" of a limb.

%create a matrix excluding time stamp data
vertex_xy = simdata(2:end);

% The vertex coordinates are stored as [x1, y1, x2, y2, .... xN, yN], but
% we'll need to separate out the x and y values. Take every other column:
vertex_x = vertex_xy(1:2:end);
vertex_y = vertex_xy(2:2:end);

%Need all 14 xy coordinate pairs for each of the separate limbs, create a
%data structure where each each field is a matrix of the xy coordiates of
%the limb.  i.e. curve_data{1} will be a (num_v_per_circ x 2) matrix with 
%the x coordiates in column 1 and the y coordinates in column 2.  Note,
%each curved section originally has 13 vertexes, but we will add the
%previous tip as the first entry

circ_x = vertex_x(1:(num_limbs * num_v_per_circ)); %x-coordinates for the curved sections
circ_y = vertex_y(1:(num_limbs * num_v_per_circ)); %y-coordinates for the curved sections


%create data structure for curved sections
for i = 1:num_limbs
    if i == 1
        circ_x_single_limb = circ_x((i-1)*num_v_per_circ+1:(i*num_v_per_circ));
        circ_y_single_limb = circ_y((i-1)*num_v_per_circ+1:(i*num_v_per_circ));
        circ_x_single_limb = [circ_x(num_limbs * num_v_per_circ),circ_x_single_limb];
        circ_y_single_limb = [circ_y(num_limbs * num_v_per_circ),circ_y_single_limb]; 
    else
        %x and y coordinates for points on the for the curved section of the ith limb
        circ_x_single_limb = circ_x((i-1)*num_v_per_circ+1:(i*num_v_per_circ));
        circ_y_single_limb = circ_y((i-1)*num_v_per_circ+1:(i*num_v_per_circ));
        circ_x_single_limb = [circ_x((i-1)*num_v_per_circ),circ_x_single_limb];
        circ_y_single_limb = [circ_y((i-1)*num_v_per_circ),circ_y_single_limb]; 

    end
    %create data structure of x and y coordinates for each limb
    curve_data{i} = [circ_x_single_limb; circ_y_single_limb];
end

stright_x = vertex_x((num_limbs * num_v_per_circ)+1:end); %x-coordinates for the curved sections
stright_y = vertex_y((num_limbs * num_v_per_circ)+1:end); %y-coordinates for the curved sections

%create data structure for stright sections
for i = 1:num_limbs
    stright_x_single_limb = stright_x((i-1)*num_v_per_flat+1:(i*num_v_per_flat));
    stright_y_single_limb = stright_y((i-1)*num_v_per_flat+1:(i*num_v_per_flat));   
    stright_x_single_limb = [curve_data{i}(1,1), stright_x_single_limb];
    stright_y_single_limb = [curve_data{i}(2,1), stright_y_single_limb]; 
    stright_data{i} = [stright_x_single_limb; stright_y_single_limb];
end

%create vectors of the x and y cooridates of the tips of the robot.  The
%tips occur as the num_v_per_circ 'th row of the curve data matricies. 
for k = 1:num_limbs
    curve_tips_x(k) = curve_data{k}(1,1);
    curve_tips_y(k) = curve_data{k}(2,1);
end
%create a matrix of tip locations col 1 = tips_x, col 2 = tips_y
curve_tips = [curve_tips_x; curve_tips_y];

%calculate the xy coordinates of the true tips of the robot
for i = 1:num_limbs
    tips_x(i) = stright_data{i}(1,num_v_per_flat+1);
    tips_y(i) = stright_data{i}(2,num_v_per_flat+1);
end
tips = [tips_x;tips_y];
end
