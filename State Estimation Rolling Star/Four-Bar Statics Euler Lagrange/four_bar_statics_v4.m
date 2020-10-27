%% RISS 2020 STATICS 4-BAR v4 %%
% Using statics, determine the necessary angles to fully define the shape
% of the robot.  This script utilizes the equler lagrange approach

%CHANGE

clc
clear all
close all

%create symbolic variables
syms kappa1 kappa2 kappa3 L1 L2 L3 a4 t1 t2 t3 m1 m2 m3 g K2 K3 g t02 t03



%% forward kineamtics
%determine bar lengths
n = 4;
K = [kappa1; kappa2; kappa3];
L = [L1; L2; L3];
% for j = 1:n-1
%     a{j} = 0;
% end

for i = 1:n-1
    a(i) = (2*sin((L(i)*K(i))/2))/K(i);
end

[T0n,Tnm1_n] = fwdkinRISS([a(1);a(2);a(3)], [0;0;0], [0;0;0], [t1;t2;t3])

%joint locations
% x0_vec = [0;0;0];
% x1_vec = [a(1)*cos(t1);
%           a(1)*sin(t1);
%           0];
% x2_vec = x1_vec + a(2)*[cos(t1+t2);sin(t1+t2);0];
% x3_vec = [-a4;0;0];
x0_vec = [0;0;0];
x1_vec = T0n{1}(1:3,4);
x2_vec = T0n{2}(1:3,4);
x3_vec = T0n{3}(1:3,4);

% Determine COM location
for i = 1:n-1
    phi = L(i)*K(i)*.5;
    A =(2*sin(phi))/(L(i)*K(i)^2);
    B = cos(phi)/K(i);
    C{i} = A-B;
end

%COM locations
P0G1 = a(1)/2*[cos(t1);sin(t1);0]+C{1}*T0n{1}(1:3,2);
P0G2 = x1_vec + a(2)/2*[cos(t1+t2);sin(t1+t2);0]+C{2}*T0n{2}(1:3,2);
P0G3 = x3_vec-a(3)/2*[cos(t1+t2+t3);sin(t1+t2+t3);0]+C{3}*T0n{3}(1:3,2);

%% Geometry
x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1));
alpha = asin((a4*sin(t1))/x);
beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
zeta = asin((x*sin(beta))/a(3));
t2_final = simplify(pi - alpha - beta);
t3_final = simplify(pi-zeta);

%% deterine Y
Y = [                          1;
    simplify(diff(t2_final,t1));
    simplify(diff(t3_final,t1))];

%% Determine gravitational forces
%determine jacobains at COM
% JVG1 = [simplify(diff(P0G1,t1)), simplify(diff(P0G1,t2)), simplify(diff(P0G1,t3))];
% JVG2 = [simplify(diff(P0G2,t1)), simplify(diff(P0G2,t2)), simplify(diff(P0G2,t3))];
% JVG3 = [simplify(diff(P0G3,t1)), simplify(diff(P0G3,t2)), simplify(diff(P0G3,t3))];

%determine g0
% g_dir = [0 -g 0].'; %gravity expressed as a 3-space vector
% g0 = -[JVG1.' JVG2.' JVG3.']*[m1.*g_dir; m2.*g_dir; m3.*g_dir];

% g_a = simplify(subs(Y.'*g0,{t2,t3},{t2_final,t3_final}))

%% Deterime spring forces
tau0 = [             0;
        K2*(pi-t2-t02);
        K3*(pi-t3-t03)];

tau_a = simplify(subs(Y.'*tau0,{t2,t3},{t2_final,t3_final}))

%% Find the residual
% R = g_a-tau_a
% dR_dt1 = simplify(diff(R,t1))








