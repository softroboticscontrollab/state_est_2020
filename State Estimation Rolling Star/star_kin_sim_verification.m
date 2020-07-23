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

%% Define constants
curve_L = 0.022; % lenghth of curved section (m)
Ls = .008; % Length of stright section
time = 0; % simulation time
num_limbs = 7; %number of limbs
num_v_per_circ = 13; % number of vertices per "circular part" of a limb
num_v_per_flat = 5; % number of vertices per "flat part" of a limb.
% ...note that total number of vertices, i.e. columns in vertex_x or vertex_y, is:
% num_limbs * (num_v_per_circ + num_v_per_flat).

%% Load data from the DER simulator
% The function to load from the simulator's CSV file is in a subdirectory
addpath( genpath('simRollingStarData') );

% Manually specify the filename. We'll assume the file is in the same
% directory as the m-file that loads it.
der_data_fname = 'simRollingStarAllNodes_2020_07_09_143539.csv';
[simdata] = load_simRollingStar_dataset(der_data_fname);

%%  Parse data from simulation at a particular time
[vertex_xy,curve_data,stright_data,circ_tips,tips] = rollingStar_time(simdata,time);

%% Do something as an example
% % The vertex coordinates are stored as [x1, y1, x2, y2, .... xN, yN], but
% we'll need to separate out the x and y values. Take every other column:
% vertex_x = vertex_xy(:, 1:2:end);
% vertex_y = vertex_xy(:, 2:2:end);

% % Plot all nodes of the robot at the first timestep.
% figure(1);
% hold on;
% plot(vertex_x(1,:), vertex_y(1,:), 'k.');
% title('All nodes at time zero');
% 
% % since the data is a bit confusing in terms of "circular" versus "flat"
% % sections, here's an example of each.
% timezero_circ_x = vertex_x(1, 1:(num_limbs * num_v_per_circ));
% timezero_circ_y = vertex_y(1, 1:(num_limbs * num_v_per_circ));
% timezero_flat_x = vertex_x(1, (num_limbs * num_v_per_circ):end);
% timezero_flat_y = vertex_y(1, (num_limbs * num_v_per_circ):end);
% 
% figure(2);
% hold on;
% plot(timezero_circ_x, timezero_circ_y, 'r.');
% title('Circular sections of the robot only');
% 
% figure(3);
% hold on;
% plot(timezero_flat_x, timezero_flat_y, 'b.');
% title('Flat sections (sticking out past the circular sections)');

%% Plot simulation data
%stright sections
% for i = 1:num_limbs
%     figure(1)
%     axis equal
%     plot(stright_data{i}(1,:),stright_data{i}(2,:),'k.')
%     hold on
% end
% %curved sections
for i = 1:num_limbs
    figure(2)
    axis equal
    plot(curve_data{i}(1,:),curve_data{i}(2,:),'rx')
    hold on
end
xlabel('x-location')
ylabel('y-location')
set(gcf,'color','w');
% %both
for i = 1:num_limbs
    figure(3)
    axis equal
    plot(curve_data{i}(1,:),curve_data{i}(2,:),'rx')
    hold on
    plot(stright_data{i}(1,:),stright_data{i}(2,:),'k.')
    hold on
end

%% Calculate curvatures of the limbs
% calculate curve fit for each of curved sections of limbs and plot
for j = 1:num_limbs
[xc,yc,Re] = circfit(curve_data{j}(1,:),curve_data{j}(2,:)); % create a curve fit for each of the limbs
K(j) = 1/Re; % record curvatures 
th = linspace(0,2*pi,100)'; % create a vector of angles to plot circle
% calcuale circle from curve fit
xe = Re*cos(th)+xc; 
ye = Re*sin(th)+yc;
% plot curve fit circles 
% figure(2) 
% hold on
% axis equal
% plot([xe;xe(1)],[ye;ye(1)])
end

% The below functions connect the points of the curved sections of 2 limbs
% to determine which of the 14 points corresponds to which limb
figure(2)
hold on 
plot(curve_data{1}(1,:),curve_data{1}(2,:))
hold on 
plot(curve_data{2}(1,:),curve_data{2}(2,:),'g')

%% dertiremine the jt angles of L-3 limbs from sim data
%(compare with statics model later on)
%calculate angle vector in form: theta = [t1p,t2,t3,t_L-3]
theta(1) = int_ang(circ_tips(:,num_limbs),circ_tips(:,1),circ_tips(:,2));
%calculate other necessary angles
for m = 2:num_limbs-3
    theta(m) = pi-int_ang(circ_tips(:,m-1),circ_tips(:,m),circ_tips(:,m+1));
end
%calcualate thetas 
theta_deg = theta*180/pi;

%% Fwd kinematics curved sections
% fwd kin function stright sections
[tip_pos,T0_n,T0_np, a] = kin_7bar(K,theta,curve_L);

%since the first leg is shifted, from 0, adjust the tip poistion
x_shift = circ_tips(1,1);
y_shift = circ_tips(2,1);
tip_pos(1,:) = tip_pos(1,:)+x_shift;
tip_pos(2,:) = tip_pos(2,:)+y_shift;

% Plot the tops of the curved sections
% figure(2)
% plot(tip_pos(1,:),tip_pos(2,:),'bx')

%% Fwd kinematics stright sections
% fwd kin function curved sections
[stright_tips,ms] = straight_kin_7bar(T0_n,T0_np,a,K,Ls);

% since the first leg is shifted, from 0, adjust the tip poistion
stright_tips(1,:) = stright_tips(1,:)+x_shift;
stright_tips(2,:) = stright_tips(2,:)+y_shift;

% plot cacluated tip locations of the stirhgt sections
for i = 1:num_limbs
    figure(3)
    axis equal
    plot(stright_tips(1,i),stright_tips(2,i),'bx')
    hold on 
end

%% calculate error