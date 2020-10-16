function [tip_pos,T0_n,T0_np, a_num] = kin_7bar(K,theta,curve_L)
%kin_7bar takes in a vector of curvature, angles and an arc length and
%calculates a matrix, tip_pos (2x7 matrix where row 1 is the x positions
%and row 2 is the tip y positions, column 1 is tip 1).  This function also
%prints the transformation matricies from fram 0 to each of the tips
%
%   Inputs: 
%
%   K = vector of curvatures of each of the 7 bars
%       in form K = [k1,k2,k3,k4,k5,k6,k1p]
%   theta = L-3 angles to fully define the robot's pose
%           in form theta = [t1p,t2,t3,t4]
%   curve_L = curve length of the curved section of the arm
%
%   Outputs:
%
%   tip_pos = the position of each of the 7 tips
%   

L = 7;  %enter the number of bars (leave as 7 for this function)

% Enter L-3 joint angles
t1 = 0;        %jt angle 1 (rad)
t1p = theta(1);  %jt angle 1' (rad)
t2 = theta(2);    %jt angle 2 (rad)
t3 = theta(3);    %jt angle 3 (rad)
t4 = theta(4);    %jt angle 4 (rad)

% Determine bar length form curvature
a_num = barcalc(K,curve_L); %determine bar lengths (m)

%% Forward kinematics with all angles
syms t6 t5 %create symbolic variables

theta = [t1,t2,t3,t4,t5,t6]; %create a vector of angles 

% Create d and alpha vectors of 0
for i=1:L-1
    alpha(i,1) = 0;
    d(i,1) = 0;
end

% Calculate kinematics from frame 0 to L-1
[T0_n,Tnm1_n] = fwdkinRISSnum(a_num(1:L-1), alpha, d, theta);

% Forward kinematics fram 0 to a1p
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
s_zeta = (x*sin(beta))/aLm1;
c_zeta = (-x^2+aLm1^2+aLm2^2)/(2*aLm1*aLm2);
zeta = atan2(s_zeta,c_zeta);
% Calculate the two unknown angles
t5_num = vpa(pi - alpha - beta - (sum(theta(1:L-3))-gamma),5);
t6_num = vpa(pi-zeta,5);

%% Find and print transformation matricies
for i = 1:L-1
    if i<L-2
%         fprintf('Transomation 0 to %2.0f\n',i)
        T0_n{i} = vpa(T0_n{i},5);
%         T0_n{i}
    else
%         fprintf('Transomation 0 to %2.0f\n',i)
        T0_n{i} = vpa(subs(T0_n{i},{t6,t5},{t6_num,t5_num}),5);
%         T0_n{i}
    end
   
end

% fprintf('Transomation 0 to 1p')
% T0_np{1}

%determine tip positions of curved sections
tip_pos = zeros(2,L);
tip_pos(1,1) = 0; %x-pos
tip_pos(2,1) = 0; %y_pos
for i = 2:L
    %x-pos
    tip_pos(1,i) = vpa(T0_n{i-1}(1,4),5);
    %y_pos
    tip_pos(2,i) = vpa(T0_n{i-1}(2,4),5);
end

end

