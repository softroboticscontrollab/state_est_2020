%%RISS 2020 STATICS 4-BAR%%

clc
clear all
close all 

%create variables
syms f2 a1 a2 a3 a1p t1 t2 t3 t1p m1 m2 m3 m1p g k1 k2 k2 k3 k1p L C m1 m2 m3 m1p

%% forward kinematics fram 0 to 2
a = [a1 a2 a3].';
alpha = [0 0 0].';
d = [0 0 0].';
theta = [0 t2 t3].';
[T0_n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta);

%% forward kinematics fram 0 to 1p
a01p = [a1p];
alpha = [0].';
d = [0].';
theta = [t1p].';
[T0_1p,T0_1] = fwdkinRISSnum(a01p, alpha, d, theta);

%% Determine Jacobians at COM
%enter necessary directions and vectors
P0_G1 = T0_n{1}(1:3,4);
P0_G2 = T0_n{2}(1:3,4);
P0_G3 = T0_n{3}(1:3,4);
P0_G1p = [0 0 0].';
z0_0 = [0 0 1].';
z0_1 = [0 0 1].';
z0_2 = [0 0 1].';

T1_3 = Tnm1_n{2}*Tnm1_n{3}
R0_1 = T0_n{1}(1:3,1:3).'
P1_G3 = R0_1*T1_3(1:3,1)

R0_2 = T0_n{2}(1:3,1:3).'
P2_G3 = R0_2*Tnm1_n{3}(1:3,1)

P1_G2 = R0_1*Tnm1_n{2}(1:3,1)

%for link 1
JV1 = simplify(cross(z0_0,P0_G1));
JV2 = [0 0 0].';
JV3 = [0 0 0].';
JV4 = [0 0 0].';

%create J matrix
JVG1 = [JV1 JV2 JV3 JV4]

%for link 2
JV1 = simplify(cross(z0_0,P0_G2));
JV2 = simplify(cross(z0_1,P1_G2));
JV3 = [0 0 0].';
JV4 = [0 0 0].';

%create J matrix
JVG2 = [JV1 JV2 JV3 JV4]

%for link 3
JV1 = simplify(cross(z0_0,P0_G3));
JV2 = simplify(cross(z0_1,P1_G3));
JV3 = simplify(cross(z0_2,P2_G3));
JV4 = [0 0 0].';

%create J matrix
JVG3 = [JV1 JV2 JV3 JV4]

%for link 1p
JV1 = (cross(z0_0,P0_G1p));
JV2 = [0 0 0].';
JV3 = [0 0 0].';
JV4 = [0 0 0].';

%create J matrix
JVG1p = [JV1 JV2 JV3 JV4]

%% create gravity vector
g_dir = [0 -g 0].';
G = -[JVG1.' JVG2.' JVG3.' JVG1p.']*[m1.*g_dir; m2.*g_dir; m3.*g_dir; m1p.*g_dir]


%% Create J2
JWG2 = [0 0 0 0;
        0 0 0 0;
        1 1 0 0;]
    
J2 = [JVG2;JWG2]

fext = [0; f2; 0; 0; 0; 0]

u = J2.'*fext