%%RISS 2020 FWD KIN 4 LEGGED ROLLING STAR%%

clc
clear all
close all 

%create variables
syms a1 a2 a3 a1p t1 t2 t3 t1p 

% x = sqrt(a1^2 + a1p^2 - 2*a1*a1p*cos(t1p-t1))
% 
% alpha = simplify(asin((a1p*sin(t1p-t1))/x))
% 
% beta = simplify(acos((x^2 + a2^2 - a3^2)/(2*x*a2)))
% 
% t2 = simplify(pi - alpha - beta)

%% forward kinematics fram 0 to 3
a = [a1 a2 a3].';
alpha = [0 0 0].';
d = [0 0 0].';
theta = [t1 t2 t3].';
[T0_n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta)

%% forward kinematics fram 0 to a1p
a = [a1p].';
alpha = [0].';
d = [0].';
theta = [t1p].';
[T0_1p,T0_1] = fwdkinRISS(a, alpha, d, theta)

P0_3 = T0_n{3}(1:3,4)
P0_1p = T0_1p{1}(1:3,4)

R3_0 = T0_n{3}(1:3,1:3).'
% 
simplify(R3_0*(P0_3-P0_1p))
% 
simplify(P0_3-P0_1p)
