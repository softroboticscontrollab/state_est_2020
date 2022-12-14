%% RISS 2020 N-BAR KINEMATICS %%
% Given the arc lengths, curvatures, material proprties and necessary 
% angles to fully define the robot's pose, determine transformation 
% matricies symbolically from the base to each of the robot's joints 
% and print them
clc
clear all
close all 

%% Sym variable names
L = 5;  %enter bar lengths
theta_max = L-3; %number of angles to fully define geometry
syms x y alpha beta zeta gamma a1 t1 a1p t1p %create symbolic variables
a(1) = a1; %initialize a vector
theta(1) = t1; %initialize theta vector

%create symbolic variables vectors of the correct length
for i=2:L-1
    bar_name = strcat('a',num2str(i));
    syms(bar_name);
    a(end+1) = bar_name;
    angle_name = strcat('t',num2str(i)); 
    syms(angle_name);
    theta(end+1) = angle_name;
end

%% Forward kinematics with all angles
%create d and alpha vectors of 0
for i=1:L-1
    alpha(i,1) = 0;
    d(i,1) = 0;
end

%calculate kinematics from frame 0 to L-1
[T0_n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta);

% forward kinematics fram 0 to a1p
[T0_np,Tnm1_np] = fwdkinRISS(a1p, 0, 0, t1p);

%% Closed loop constraints
xLm3 = T0_n{L-3}(1,4);
yLm3 = T0_n{L-3}(2,4);
aLm1 = a(L-1);
aLm2 = a(L-2);
y = sqrt(xLm3^2+yLm3^2);
gamma = atan2(yLm3,xLm3);
x = sqrt(y^2 + a1p^2 - 2*y*a1p*cos(t1p-gamma));
alpha = asin((a1p*sin(t1p-gamma))/x);
beta = acos((x^2 + aLm2^2 - aLm1^2)/(2*x*aLm2));
zeta = asin((x*sin(beta))/aLm1);
% Calculate the two unknown angles
tLm2 = pi - alpha - beta - sum(theta(1:L-3));
tLm1 = pi-zeta;

% print symbolic solution
for i = 1:L-1
    if i<L-2
        fprintf('Transomation 0 to %2.0f\n',i)
        T0_n{i}
    else
        fprintf('Transomation 0 to %2.0f\n',i)
        (subs(T0_n{i},{theta(L-1),theta(L-2)},{tLm2,tLm1}))
    end
 end
 
fprintf('Transomation 0 to 1p')
T0_np{1}