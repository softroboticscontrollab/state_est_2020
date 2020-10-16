%% RISS 2020 STATICS 4-BAR v4 %%
% Using statics, determine the necessary angles to fully define the shape
% of the robot.  This script utilizes the equler lagrange approach

%CHANGE

clc
clear all
close all 

%create symbolic variables
syms a1 a2 a3 a4 t1 t2 t3 m1 m2 m3 g k2 k3 g t0

%% forward kineamtics
%joint locations
x0_vec = [0;0;0];
x1_vec = [a1*cos(t1);
          a1*sin(t1);
          0];
x2_vec = x1_vec + a2*[cos(t1+t2);sin(t1+t2);0];
x3_vec = [-a4;0;0];

%COM locations  
P0G1 = a1/2*[cos(t1);sin(t1);0];
P0G2 = x1_vec + a2/2*[cos(t1+t2);sin(t1+t2);0];
P0G3 = x3_vec-a3/2*[cos(t1+t2+t3);sin(t1+t2+t3);0];

%% Geometry
x = sqrt(a1^2 + a4^2 - 2*a1*a4*cos(pi-t1));
alpha = asin((a4*sin(t1))/x);
beta = acos((x^2 + a2^2 - a3^2)/(2*x*a2));
zeta = asin((x*sin(beta))/a3);
t2_final = simplify(pi - alpha - beta);
t3_final = simplify(pi-zeta);

%% deterine Y
Y = [                          1;
     simplify(diff(t2_final,t1));
     simplify(diff(t3_final,t1))];

%% Determine gravitational forces
%determine jacobains at COM
JVG1 = [simplify(diff(P0G1,t1)), simplify(diff(P0G1,t2)), simplify(diff(P0G1,t3))];
JVG2 = [simplify(diff(P0G2,t1)), simplify(diff(P0G2,t2)), simplify(diff(P0G2,t3))];
JVG3 = [simplify(diff(P0G3,t1)), simplify(diff(P0G3,t2)), simplify(diff(P0G3,t3))];

%determine g0
g_dir = [0 -g 0].'; %gravity expressed as a 3-space vector
g0 = -[JVG1.' JVG2.' JVG3.']*[m1.*g_dir; m2.*g_dir; m3.*g_dir]; 

g_a = simplify(subs(Y.'*g0,{t2,t3},{t2_final,t3_final}))

%% Deterime spring forces
tau0 = [            0;
        k2*(pi-t2-t0);
        k3*(pi-t3-t0)];
    
tau_a = simplify(subs(Y.'*tau0,{t2,t3},{t2_final,t3_final}))

%% Find the residual
R = g_a-tau_a
dR_dt1 = simplify(diff(R,t1))







     