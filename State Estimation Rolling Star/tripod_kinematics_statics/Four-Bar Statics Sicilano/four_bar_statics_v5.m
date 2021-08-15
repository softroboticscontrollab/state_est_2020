%% RISS 2020 STATICS 4-BAR v5 %%
% Using statics, determine the necessary angles to fully define the shape
% of the robot.  This script utilizes the equler lagrange approach

%% Prep the workspace
clc
clear all
close all

%% Create symbolic variables
syms kappa1 kappa2 kappa3 L1 L2 L3 a1 a2 a3 a4 t1 t2 t3 m1 m2 m3 g k2 k3 g t02 t03 F

%% Determine setup configuration
COM_on_bar = true;
gravity = false;
force_to_right = false;

%% Forward kineamtics
% Determine bar lengths
n = 4;
if (COM_on_bar)
    a = [a1,a2,a3];
    % Run forward kineamtics algorithm
    [T0n,Tnm1_n] = fwdkinRISS([a(1);a(2);a(3)], [0;0;0], [0;0;0], [t1;t2;t3]);
    % COM locations
    P0G1 = a(1)/2*[cos(t1);sin(t1);0];
    P0G2 = T0n{1}(1:3,4) + a(2)/2*[cos(t1+t2);sin(t1+t2);0];
    P0G3 = T0n{3}(1:3,4)-a(3)/2*[cos(t1+t2+t3);sin(t1+t2+t3);0];
else
    K = [kappa1; kappa2; kappa3];
    L = [L1; L2; L3];
    for i = 1:n-1
        a(i) = (2*sin((L(i)*K(i))/2))/K(i);
    end
    % Run forward kineamtics algorithm
    [T0n,Tnm1_n] = fwdkinRISS([a(1);a(2);a(3)], [0;0;0], [0;0;0], [t1;t2;t3]);
    % Determine COM location
    for i = 1:n-1
        phi = L(i)*K(i)*.5;
        A =(2*sin(phi))/(L(i)*K(i)^2);
        B = cos(phi)/K(i);
        C{i} = A-B;
    end
    % COM locations
    P0G1 = a(1)/2*[cos(t1);sin(t1);0];
    P0G2 = T0n{1}(1:3,4); + a(2)/2*[cos(t1+t2);sin(t1+t2);0];
    P0G3 = T0n{3}(1:3,4)-a(3)/2*[cos(t1+t2+t3);sin(t1+t2+t3);0];
end

%% Geometry
x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1));
alpha = asin((a4*sin(t1))/x);
beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
% zeta = asin((x*sin(beta))/a(3));
t2_final = (pi - alpha - beta);
t3_final = (pi-zeta);

%% deterine Y
Y = [                  1;
    (diff(t2_final,t1));
    (diff(t3_final,t1))]

%% Determine gravitational forces
%determine jacobains at COM
JVG1 = [(diff(P0G1,t1)), (diff(P0G1,t2)), (diff(P0G1,t3))];
JVG2 = [(diff(P0G2,t1)), (diff(P0G2,t2)), (diff(P0G2,t3))];
JVG3 = [(diff(P0G3,t1)), (diff(P0G3,t2)), (diff(P0G3,t3))];

% determine g0 and external forces (putting external force in go since it is
% also acting on COM
g_dir = [0 -g 0].'; %gravity expressed as a 3-space vector
if (gravity==true) &&  (force_to_right == true)
    F3 = [F -g 0].';
    g0 = -[JVG1.' JVG2.' JVG3.']*[m1.*g_dir; m2.*g_dir; F3];
elseif (gravity==true) &&  (force_to_right == false)
    g_dir = [0 -g 0].'; %gravity expressed as a 3-space vector
    g0 = -[JVG1.' JVG2.' JVG3.']*[m1.*g_dir; m2.*g_dir; m3.*g_dir];
    
elseif (gravity==false) &&  (force_to_right == true)
    zero_vec = [0 0 0].';
    g0 = -[JVG1.' JVG2.' JVG3.']*[zero_vec; zero_vec; [F 0 0].'];
else
    g0 = [0 0 0].';
end
g_a = (subs(Y.'*g0,{t2,t3},{t2_final,t3_final}))

%% Deterime spring forces
tau0 = [             0;
    k2*(pi-t2-t02);
    k3*(pi-t3-t03)]

tau_a = (subs(Y.'*tau0,{t2,t3},{t2_final,t3_final}))

%% Find the residual
R = g_a-tau_a
dR_dt1 = (diff(R,t1))








