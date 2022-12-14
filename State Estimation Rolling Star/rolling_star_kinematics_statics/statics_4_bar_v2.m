%%RISS 2020 STATICS 4-BAR v2%%

clc
clear all
close all 

%create variables
syms f2y a1 a2 a3 a1p t0 t1 t2 t3 t1p m1 m2 m3 m1p g k1 k2 k2 k3 k1p 

%% forward kinematics fram 0 to 3
a = [a1 a2 a3].';
alpha = [0 0 0].';
d = [0 0 0].';
theta = [t1 t2 t3].';
[T0_n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta);

%% forward kinematics fram 0 to 1p
[T0_np,Tnm1_np] = fwdkinRISS(a1p, 0, 0, t1p);

%% Determine Jacobians at COM
%enter necessary directions and vectors
z00 = T0_n{1}(1:3,1);
P01 = T0_n{1}(1:3,4);
P01p = T0_np{1}(1:3,4);
P02 = T0_n{2}(1:3,4)
x20 = T0_n{2}(1:3,1);
x30 = T0_n{3}(1:3,1);
x1p0 = T0_np{1}(1:3,1);
P0G1 = (a1/2).*[cos(t1);sin(t1);0];
P0G2 = P01+(a2/2)*x20;
% P0G3 = P01p - (a3/2)*x30;
P0G3 = P02 + (a3/2)*x30
P0G1p = (a1p/2)*x1p0;

%for link 1
JV1 = simplify(diff(P0G1,t1));
JV2 = simplify(diff(P0G1,t2));
JV3 = simplify(diff(P0G1,t3));
JV4 = simplify(diff(P0G1,t1p));

%create J matrix
JVG1 = [JV1 JV2 JV3 JV4];

%for link 2
JV1 = simplify(diff(P0G2,t1));
JV2 = simplify(diff(P0G2,t2));
JV3 = simplify(diff(P0G2,t3));
JV4 = simplify(diff(P0G2,t1p));

%create J matrix
JVG2 = [JV1 JV2 JV3 JV4];

%for link 3
JV1 = simplify(diff(P0G3,t1));
JV2 = simplify(diff(P0G3,t2));
JV3 = simplify(diff(P0G3,t3));
JV4 = simplify(diff(P0G3,t1p));

%create J matrix
JVG3 = [JV1 JV2 JV3 JV4];

%for link 1p
JV1 = simplify(diff(P0G1p,t1));
JV2 = simplify(diff(P0G1p,t2));
JV3 = simplify(diff(P0G1p,t3));
JV4 = simplify(diff(P0G1p,t1p));

%create J matrix
JVG1p = [JV1 JV2 JV3 JV4];

%% create gravity vector
g_dir = [0 -g 0].';
G = -[JVG1.' JVG2.' JVG3.' JVG1p.']*[m1.*g_dir; m2.*g_dir; m3.*g_dir; m1p.*g_dir]

%% external force analysis
z_vec = [0;0;0];
Jtop = [simplify(diff(P0G1,t1)) z_vec z_vec z_vec];
Jbottom = [z00 z_vec z_vec z_vec];
J1 = [Jtop;Jbottom];
gamma = [0;f2y;0;0;0;0];
statics_ext = J1.'*gamma

%% geometry
x = sqrt(a1^2 + a1p^2 - 2*a1*a1p*cos(t1p-t1));
alpha = asin((a1p*sin(t1p-t1))/x);
beta = acos((x^2 + a2^2 - a3^2)/(2*x*a2));
t2_final = pi - alpha - beta;
s_zeta = (x*sin(beta))/a3;
c_zeta = (-x^2+a3^2+a2^2)/(2*a3*a2);
zeta = atan2(s_zeta,c_zeta);
t3_final = pi-zeta;

%% spring torqes
tau1 = k1*(t0 - t1');
tau2 = k2*(t2_final-t0);
tau3 = k3*(t3_final-t0);
tau1p = k1p*(t0-(2*pi-t1p-alpha-beta-zeta));
tau = [t1,t2,t3,t1p]

%% solution in terms of t1p
final = subs(G+statics_ext,{t1,t2,t3},{0,t2_final,t3_final});
final = (subs(final,{t1,t2,t3},{0,t2_final,t3_final}))

final_nosub = simplify(G+statics_ext)
