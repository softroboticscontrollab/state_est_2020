%%RISS 2020 N-BAR KINEMATICS%%

clc
clear all
close all 

%% sym variable names
L = 7;  %enter bar lengths
theta_max = L-3; %number of angles to fully define geometry
syms x y alpha beta zeta gamma a1 t1 a1p t1p %create symbolic variables
a(1) = a1; %initialize a vector
theta(1) = t1; %initialize theta vector

%create symbolic variables
for i=2:L-1
    bar_name = strcat('a',num2str(i));
    syms(bar_name);
    a(end+1) = bar_name;
    angle_name = strcat('t',num2str(i)); 
    syms(angle_name);
    theta(end+1) = angle_name;
end

%% forward kinematics with all angles
%create d and alpha vectors of 0
for i=1:L-1
    alpha(i,1) = 0;
    d(i,1) = 0;
end

%calculate kinematics
[T0_n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta);


%% forward kinematics fram 0 to a1p
% a = [a1p].';
% alpha = [0].';
% d = [0].';
% theta = [t1p].';

%calculate kinematics
[T0_np,Tnm1_np] = fwdkinRISS(a1p, 0, 0, t1p);

%% closed loop constraints
xLm2 = T0_n{L-2}(1,4);
yLm2 = T0_n{L-2}(2,4);
aLm1 = a(L-1);
aLm2 = a(L-2);
y = sqrt(xLm2^2+yLm2^2);
gamma = atan2(yLm2,xLm2);

x = sqrt(y^2 + a1p^2 - 2*y*a1p*cos(t1p-gamma));
alpha = asin((a1p*sin(t1p-gamma))/x);
beta = acos((x^2 + aLm2^2 - aLm1^2)/(2*x*a2));

tLm2 = pi - alpha - beta - sum(theta,[1,L-2]);
zeta = asin((x*sin(beta))/aLm1);
tLm1 = pi-zeta;

%% print symbolic solution
for i = 1:L-1
    if i<L-2
        fprintf('Transomation %2.0f\n',i)
        T0_n{i}
    else
        fprintf('Transomation %2.0f\n',i)
        simplify(subs(T0_n{i},{theta(L-1),theta(L-2)},{tLm2,tLm1}))
    end
   

    end
%     
% T0_n{1}
% 
% subs(theta(L-2),tLm2)
% subs(theta(L-1),tLm1)



% angle_name = strcat('t',num2str(i)); 
% syms(angle_name);
% theta(end+1) = angle_name;
% 
% t = pi - alpha - beta;


