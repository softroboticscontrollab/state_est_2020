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
time = 0.24; % simulation time
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


simtime = simdata(:,1);

for i = 1:length(simtime)
    sim_xy{i} = zeros(2,126);
    sim_COM_xy{i} = zeros(2,1);
    model_COM_xy{i} = zeros(2,1);
    curve_xy{i} = zeros(2,7);
    straight_xy{i} = zeros(2,7);
end

for q = 1:length(simtime)
    time = simtime(q);
    
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
e_center(q) = sqrt((simCOM(1)-modelCOM(1))^2+(simCOM(2)-modelCOM(2))^2);
e_curve_tips = sqrt((circ_tips(1,:)-tip_pos(1,:)).^2+(circ_tips(2,:)-tip_pos(2,:)).^2);
e_stright_tips = sqrt((tips(1,:)-straight_tips(1,:)).^2+(tips(2,:)-straight_tips(2,:)).^2);

%calculate SEE
SEE_curve(q) = sum(e_curve_tips.^2)/(length(e_curve_tips)-2);
SEE_straight(q) = sum(e_stright_tips.^2)/(length(e_stright_tips)-2);

% % alternate SEE, avg
AvgE_curve(q) = sum(e_curve_tips)/length(e_curve_tips);
AvgE_straight(q) = sum(e_stright_tips)/length(e_stright_tips);

%% Video data
sim_xy{q}(1,:) = vertex_x;
sim_xy{q}(2,:) = vertex_y;

sim_COM_xy{q}(1,:) = simCOM(1);
sim_COM_xy{q}(2,:) = simCOM(2);

model_COM_xy{q}(1,:) = modelCOM(1);
model_COM_xy{q}(2,:) = modelCOM(2);

curve_xy{q}(1,:) = tip_pos(1,:);
curve_xy{q}(2,:) = tip_pos(2,:);

straight_xy{q}(1,:) = straight_tips(1,:);
straight_xy{q}(2,:) = straight_tips(2,:);


end


%calculate % error
h = .05058;
curve_percent = 100*AvgE_curve/h;
straight_percent = 100*AvgE_straight/h;
COM_percent = 100*e_center/h;

SEE_center = sum(e_center.^2)/(length(e_center)-2); %must do after iteration

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


% % % % USE THIS FOR PAPER
% for i = 1:num_limbs
%     figure(3)
%     axis equal
%     curve_sim_plot = plot(curve_data{i}(1,:),curve_data{i}(2,:),'k.')
%     hold on
%     straight_sim_plot = plot(straight_data{i}(1,:),straight_data{i}(2,:),'k.')
%     hold on
% end
% % % % Plot model data - curve tips % % %
% figure(3)
% curve_tips_plot = plot(tip_pos(1,:),tip_pos(2,:),'x','MarkerSize',8,'color',[.2 .6 .9],'LineWidth',1.5)
%  
% % % % Plot model data - straight tips % % %
% for i = 1:num_limbs
%     figure(3)
%     axis equal
%     stright_tips_plot = plot(straight_tips(1,i),straight_tips(2,i),'rx','MarkerSize',8,'LineWidth',1.5)
%     hold on
% end
%  
% % % % Plot complete COM from simulation % % %
% figure(3)
% COM_sim_plot = plot(simCOM(1),simCOM(2),'ko','MarkerSize',8,'LineWidth',1.5)
%  
% %% % % Plot complete COM from model % % %
% figure(3)
% COM_model_plot = plot(modelCOM(1),modelCOM(2),'o','MarkerSize',8,'color',[0 .7 0],'LineWidth',1.5)
% legend([curve_sim_plot,COM_sim_plot,COM_model_plot,stright_tips_plot,curve_tips_plot],'Sim Data','COM Sim','COM Model','Straight Tips','Curve Tips','Location','southeast')
% xlabel('x-position (m)')
% ylabel('y-position (m)')
% axis([-.03 .033 -.002 .052])
% % % % END



% The below functions connect the points of the curved sections of 2 limbs
% to determine which of the 14 points corresponds to which limb
% figure(2)
% hold on
% plot(curve_data{1}(1,:),curve_data{1}(2,:))
% hold on
% plot(curve_data{2}(1,:),curve_data{2}(2,:),'g')

% % % Plot model data - curve tips % % %
% figure(3)
% curve_tips_plot = plot(tip_pos(1,:),tip_pos(2,:),'x','MarkerSize',8,'color',[.2 .6 .9],'LineWidth',1.5)

% % % Plot model data - straight tips % % %
% for i = 1:num_limbs
%     figure(3)
%     axis equal
%     stright_tips_plot = plot(straight_tips(1,i),straight_tips(2,i),'rx','MarkerSize',8,'LineWidth',1.5)
%     hold on
% end
%legend([curve_sim_plot,straight_sim_plot,stright_tips_plot],'Curve data','Straight data','Straight tips')

% % % Plot complete COM from simulation % % %
% figure(3)
% COM_sim_plot = plot(simCOM(1),simCOM(2),'ko','MarkerSize',8,'LineWidth',1.5)

% % % Plot complete COM from model % % %
% figure(3)
% COM_model_plot = plot(modelCOM(1),modelCOM(2),'mo','MarkerSize',8,'LineWidth',1.5)
% legend([curve_sim_plot,COM_sim_plot,COM_model_plot,stright_tips_plot,curve_tips_plot],'Sim Data','COM Sim','COM Model','Straight Tips','Curve Tips','Location','southeast')
% xlabel('x-position (m)')
% ylabel('y-position (m)')
% axis([-.03 .033 -.002 .052])

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
% plot(simtime,COM_percent)
% xlabel('Time(s)')
% ylabel('Error/h (%)')
% axis([0 .62 0 1.8])

% % % Plot SEE curve tips of robot over all sim times % % %
% figure(5)
% plot(simtime,curve_percent)
% xlabel('Time(s)')
% ylabel('(Avg Error) / h (%)')
% axis([0 .62 0 1.8])

% % % Plot SEE stright tips of robot over all sim times % % %
% figure(6)
% plot(simtime,straight_percent)
% xlabel('Time(s)')
% ylabel('(Avg Error) / h (%)')
% axis([0 .62 0 1.8])

% % % Plot all error on one % % % USE THIS
figure(7)
together = plot(simtime,COM_percent,'color',[0 .7 0],'LineWidth',1.5)
hold on 
plot(simtime,straight_percent,'r','LineWidth',1.5)
hold on
plot(simtime,curve_percent,'color',[.2 .6 .9],'LineWidth',1.5)
xlabel('Time(s)')
ylabel('Normalized Error (%)')
legend('Error_{com}','Error_{tips}','Error_{curve}')
axis([0 .625 0 1.8])
% % % STOP HERE



% % % Create Error Subplot % % % 
% figure
% plot(simtime,SEE_curve)
% xlabel('Time(s)')
% ylabel('Tip Location SEE (m)')

% figure
% % subplot(2,1,2)
% plot(simtime,e_center)
% xlabel('Time(s)')
% ylabel('COM Error (m)')

 
% figure
% % subplot(2,1,2)
% plot(simtime,SEE_straight)
% xlabel('Time(s)')
% ylabel('Tip Location SEE (m)')

%% Create Movie
% for i = 1:length(simtime)
%     figure(1)
%     vertexes = plot(sim_xy{i}(1,:),sim_xy{i}(2,:),'k.');
%     hold on
%     simCOM =  plot(sim_COM_xy{i}(1,:),sim_COM_xy{i}(2,:),'ko');
%     hold on
%     modelCOM = plot(model_COM_xy{i}(1,:),model_COM_xy{i}(2,:),'o','MarkerSize',8,'color',[0 .7 0],'LineWidth',1.5);
%     hold on
%     straightTips = plot(straight_xy{i}(1,:),straight_xy{i}(2,:),'rx','MarkerSize',8,'LineWidth',1.5);
%     hold on
%     curveTips = plot(curve_xy{i}(1,:),curve_xy{i}(2,:),'x','MarkerSize',8,'color',[.2 .6 .9],'LineWidth',1.5);
%     axis([-.032 .033 -.005 .052])
% %     axis equal
%     xlabel('x-position (m)')
%     ylabel('y-position (m)')
%     time = sprintf('%.3f',simtime(i))
%     title(['Sim Time: ',time,' s'])
%     legend('Sim Data','COM Sim','COM Model','Straight Tips','Curve Tips','Location','southeast')
%     set(gcf,'color','w');
%     drawnow
%     F(i) = getframe(gcf);
%     pause(.001)
%     delete(vertexes)
%     delete(simCOM)
%     delete(modelCOM)
%     delete(straightTips)
%     delete(curveTips)
% end
% 
% video = VideoWriter('Sim.mp4','MPEG-4');
% video.FrameRate = 20;
% open(video)
% writeVideo(video,F);
% close(video)

%COM - model










