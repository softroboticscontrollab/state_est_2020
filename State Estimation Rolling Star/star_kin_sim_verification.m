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
time = .24; % simulation time
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


% simtime = simdata(:,1);
% 
% for q = 1:length(simtime)
%     time = simtime(q);
    
%%  Parse data from simulation at a particular time
[vertex_xy,curve_data,straight_data,circ_tips,tips] = rollingStar_time(simdata,time);

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

%curved sections
% for i = 1:num_limbs
%     figure(2)
%     axis equal
%     plot(curve_data{i}(1,:),curve_data{i}(2,:),'rx')
%     hold on
% end
% xlabel('x-location')
% ylabel('y-location')
% set(gcf,'color','w');
% 
% %both
% for i = 1:num_limbs
%     figure(3)
%     axis equal
%     plot(curve_data{i}(1,:),curve_data{i}(2,:),'rx')
%     hold on
%     plot(straight_data{i}(1,:),straight_data{i}(2,:),'k.')
%     hold on
% end

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

%% Fwd kinematics stright sections
% fwd kin function curved sections
[straight_tips,ms,Pc] = straight_kin_7bar(T0_n,T0_np,a,K,Ls);

% since the first leg is shifted, from 0, adjust the tip poistion
straight_tips(1,:) = straight_tips(1,:)+x_shift;
straight_tips(2,:) = straight_tips(2,:)+y_shift;

%% Calculate COM
% Simulation - COM of curved sections
COM_curve = zeros(2,num_limbs);
for i = 1:num_limbs-1
    x = a(i)/(curve_L*K(i));
    P_COM = Pc{i}+x*T0_n{i}(1:3,2);
    COM_curve(1,i) = P_COM(1)+x_shift;
    COM_curve(2,i) = P_COM(2)+y_shift;
end
x = a(num_limbs)/(curve_L*K(num_limbs));
P_COM = Pc{num_limbs}-x*T0_np{1}(1:3,2);
COM_curve(1,num_limbs) = P_COM(1)+x_shift;
COM_curve(2,num_limbs) = P_COM(2)+y_shift;

%Simulation - COM of stright sections
%stright tips start one limb to the right of curve tips...
COM_straight = zeros(2,num_limbs);
COM_straight(1,1) = (straight_tips(1,num_limbs)+tip_pos(1,1))/2;
COM_straight(2,1) = (straight_tips(2,num_limbs)+tip_pos(2,1))/2;
for i = 2:num_limbs
    COM_straight(1,i) = (straight_tips(1,i-1)+tip_pos(1,i))/2;
    COM_straight(2,i) = (straight_tips(2,i-1)+tip_pos(2,i))/2;
end

%calculate total COM for model
modelCOM = zeros(2,1);
modelCOM(1) = (sum(COM_curve(1,:))+sum(COM_straight(1,:)))/(length(COM_curve(1,:))+length(COM_straight(1,:)));
modelCOM(2) = (sum(COM_curve(2,:))+sum(COM_straight(2,:)))/(length(COM_curve(2,:))+length(COM_straight(2,:)));

% Simulation
vertex_x = vertex_xy(1:2:end);
vertex_y = vertex_xy(2:2:end);
simCOM = zeros(2,1);
simCOM(1) = sum(vertex_x)/length(vertex_x);
simCOM(2) = sum(vertex_y)/length(vertex_y);

%% Calculate error
% Note:
% simCOM = COM of entire robot from simulation data
% modelCOM = COM of entire robot from model
% circ_tips = location of curve tips from simulation data
% tip_pos = location of curve tips from model
% tips = location of stright tips from simulation data
% staight_tips = location of stright tips from model

% calculate error
% e_center(q) = sqrt((simCOM(1)-modelCOM(1))^2+(simCOM(2)-modelCOM(2))^2);
% e_curve_tips = sqrt((circ_tips(1,:)-tip_pos(1,:)).^2+(circ_tips(2,:)-tip_pos(2,:)).^2);
% e_stright_tips = sqrt((tips(1,:)-straight_tips(1,:)).^2+(tips(2,:)-straight_tips(2,:)).^2);

%calculate SEE
% SEE_curve(q) = sum(e_curve_tips.^2)/(length(e_curve_tips)-2);
% SEE_straight(q) = sum(e_stright_tips.^2)/(length(e_stright_tips)-2);

% end
% SEE_center = sum(e_center.^2)/(length(e_center)-2); %must do after iteration

%% Create plots
% % % Plot simulation data % % %
%stright sections
% for i = 1:num_limbs
%     figure(1)
%     axis equal
%     plot(stright_data{i}(1,:),stright_data{i}(2,:),'k.')
%     hold on
% end
%curved sections
% for i = 1:num_limbs
%     figure(2)
%     axis equal
%     plot(curve_data{i}(1,:),curve_data{i}(2,:),'rx')
%     hold on
% end
% xlabel('x-location')
% ylabel('y-location')
% set(gcf,'color','w');
% both
for i = 1:num_limbs
    figure(3)
    axis equal
    plot(curve_data{i}(1,:),curve_data{i}(2,:),'rx')
    hold on
    plot(straight_data{i}(1,:),straight_data{i}(2,:),'k.')
    hold on
end
% The below functions connect the points of the curved sections of 2 limbs
% to determine which of the 14 points corresponds to which limb
% figure(2)
% hold on
% plot(curve_data{1}(1,:),curve_data{1}(2,:))
% hold on
% plot(curve_data{2}(1,:),curve_data{2}(2,:),'g')

% % % Plot model data - curve tips % % %
% figure(2)
% plot(tip_pos(1,:),tip_pos(2,:),'bx')

% % % Plot model data - straight tips % % %
for i = 1:num_limbs
    figure(3)
    axis equal
    plot(straight_tips(1,i),straight_tips(2,i),'bx')
    hold on
end

% % % Plot complete COM from model % % %
% figure(3)
% plot(modelCOM(1),modelCOM(2),'kx')

% % % Plot complete COM from simulation % % %
% figure(3)
% plot(simCOM(1),simCOM(2),'bx')

% % % Plot COM of curve sections from model % % %
% for i = 1:num_limbs
%     figure(3)
%     axis equal
%     plot(COM_curve(1,i),COM_curve(2,i),'gx')
%     hold on
% end

% % % Plot COM of straight sections from model % % %
% for i = 1:num_limbs
%     figure(3)
%     axis equal
%     plot(COM_straight(1,i),COM_straight(2,i),'cx')
%     hold on
% end

% % % Plot error COM of robot over all sim times % % %
% figure(4)
% plot(simtime,e_center)
% title('Error COM of robot')
% xlabel('time(s)')
% ylabel('error')

% % % Plot SEE curve tips of robot over all sim times % % %
% figure(5)
% plot(simtime,SEE_curve)
% title('Error SEE curve tips')
% xlabel('time(s)')
% ylabel('SEE')

% % % Plot SEE stright tips of robot over all sim times % % %
% figure(6)
% plot(simtime,SEE_straight)
% title('Error SEE stright tips')
% xlabel('time(s)')
% ylabel('SEE')