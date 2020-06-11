%%RISS 2020 FWD KIN 5 LEGGED ROLLING STAR%%

clc
clear all
close all 

%create variables
syms a1 a2 a3 a1p a2p t1 t2 t3 t1p t2p

%% forward kinematics fram 0 to 3
a = [a1 a2 a3].';
alpha = [0 0 0].';
d = [0 0 0].';
theta = [0 t2 t3].';
[T0_n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta)

%% forward kinematics fram 0 to a2p
a = [a1p a2p].';
alpha = [0 0].';
d = [0 0].';
theta = [t1p t2p].';
[T0_np,Tnm1_np] = fwdkinRISS(a, alpha, d, theta)

P0_3 = T0_n{3}(1:3,4)
P0_2p = T0_np{2}(1:3,4)

R3_0 = T0_n{3}(1:3,1:3).'

simplify(R3_0*(P0_3-P0_2p))

z = simplify(P0_3-P0_2p)


simplify(expand(((a3*cos(t2 + t3) - a2p*cos(t1p + t2p) + a2*cos(t2) - a1p*cos(t1p))-(a3*sin(t2 + t3) - a2p*sin(t1p + t2p) + a2*sin(t2) - a1p*sin(t1p)))^2))