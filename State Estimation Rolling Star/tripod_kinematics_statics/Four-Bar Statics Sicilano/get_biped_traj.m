function [kappa, q, a, s, masses, times] = get_biped_traj(der_data_fname, t_0, max_t)
%get_biped_traj Read DER csv data from the biped robot, calculate the
%reduced state, for a whole time trajectory
%
%   Inputs:
%       der_data_fname: filename to pass to the CSV parser
%       t0, max_t: the time frame from within the file to read
%
%   Outputs:
%       kappa, q, a: time series of curvatures, angles, and bar lengths
%       s, masses: 3x1 vectors of arc lengths and limb masses
%       times: Tx1 vector of all the timepoints we read

simdata = csvread(der_data_fname, 5, 0);
% We want all the timepoints. The function takes in time in sec, so we need
% to specify a dt to span out everything we want
dt = 0.005;
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
% n_per_limb = 13;% number of nodes for one limb

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

end

