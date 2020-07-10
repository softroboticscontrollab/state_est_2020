%% RISS 2020 7-BAR KINEMATICS NUMERIC SOLVER %%
% Given the arc lengths, curvatures, material proprties and necessary 
% angles to fully define the robot's pose, determine transformation 
% matricies from the base to each of the robot's joints

clc
clear all
close all 

%% Enter physical paramteres and find bar lengths
L = 7;  %enter the number of bars (leave as 7 for this script)

% Enter L-3 joint angles
t1 = 0 * pi/180;        %jt angle 1 (rad)
t1p = 128.57 * pi/180;  %jt angle 1' (rad)
t2 = 51.43 * pi/180;    %jt angle 2 (rad)
t3 = 51.43 * pi/180;    %jt angle 3 (rad)
t4 = 51.43 * pi/180;    %jt angle 4 (rad)

% Determine bar length form curvature
curve_L = .05;  %length of curved arm (m)
k1 = 15;        %curvature of arm 1 (m^-1)
k2 = 15;        %curvature of arm 2 (m^-1)
k3 = 15;        %curvature of arm 3 (m^-1)
k4 = 15;        %curvature of arm 4 (m^-1)
k5 = 15;        %curvature of arm 5 (m^-1)
k6 = 15;        %curvature of arm 6 (m^-1)
k1p = 15;       %curvature of arm 1p (m^-1)

K = [k1,k2,k3,k4,k5,k6,k1p]; %create vector of curvatures
a_num = barcalc(K,curve_L); %determine bar lengths (m)

%% Forward kinematics with all angles
syms x y alpha beta zeta gamma t6 t5 %create symbolic variables
theta = [t1,t2,t3,t4,t5,t6]; %create a vector of angles 

%create d and alpha vectors of 0
for i=1:L-1
    alpha(i,1) = 0;
    d(i,1) = 0;
end

%calculate kinematics
[T0_n,Tnm1_n] = fwdkinRISSnum(a_num(1:L-1), alpha, d, theta);

%% Forward kinematics fram 0 to a1p
[T0_np,Tnm1_np] = fwdkinRISSnum(a_num(L), 0, 0, t1p);

%% Closed loop constraints
xLm3 = vpa(T0_n{L-3}(1,4),5);
yLm3 = vpa(T0_n{L-3}(2,4),5);
aLm1 = a_num(L-1);
aLm2 = a_num(L-2);
a1p = a_num(7);
y = vpa(sqrt(xLm3^2+yLm3^2),5);
gamma = atan2(yLm3,xLm3);
x = vpa(sqrt(y^2 + a1p^2 - 2*y*a1p*cos(t1p-gamma)),5);
alpha = vpa(asin((a1p*sin(t1p-gamma))/x),5);
beta = acos((x^2 + aLm2^2 - aLm1^2)/(2*x*aLm2));
t5_num = vpa(pi - alpha - beta - (sum(theta(1:L-3))-gamma),5);
s_zeta = (x*sin(beta))/aLm1;
c_zeta = (-x^2+aLm1^2+aLm2^2)/(2*aLm1*aLm2);
zeta = atan2(s_zeta,c_zeta);
t6_num = vpa(pi-zeta,5);

%% Find and print solution
for i = 1:L-1
    if i<L-2
        fprintf('Transomation 0 to %2.0f\n',i)
        T0_n{i} = vpa(T0_n{i},5);
        T0_n{i}
    else
        fprintf('Transomation 0 to %2.0f\n',i)
        T0_n{i} = vpa(subs(T0_n{i},{t6,t5},{t6_num,t5_num}),5);
        T0_n{i}
    end
   
end

fprintf('Transomation 0 to 1p')
T0_np{1}

% theta_num_deg = vpa(subs(theta,{t5,t6},{t5_num,t6_num})*180/pi,6)
% theta_num_rad = vpa(subs(theta,{t5,t6},{t5_num,t6_num}),6)

%% Plot solution
x(1) = 0;
y(1) = 0;

for i = 1:length(T0_n)
    x(i+1) = T0_n{i}(1,4);
    y(i+1) = T0_n{i}(2,4); 
end

x(i+2) = 0;
y(i+2) = 0;

p1 = plot(x(1:2),y(1:2),'g-',x(2:3),y(2:3),'b-',x(3:4),y(3:4),'-',x(4:5),y(4:5),'r-',x(5:6),y(5:6),'k-',x(6:7),y(6:7),'m-',x(7:8),y(7:8),'c-');
width = 3;

for i = 1:length(p1)
    set(p1(i),{'LineWidth'},{width})
end

legend('arm 1','arm 2','arm 3','arm 4','arm 5','arm 6','arm 1p')
axis equal
xlabel('x-location (m)')
ylabel('y-location (m)')
set(gcf,'color','w');
