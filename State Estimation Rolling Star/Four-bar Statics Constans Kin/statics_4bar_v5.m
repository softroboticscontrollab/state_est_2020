%% RISS 2020 STATICS 4-BAR v5 %%
% Using statics, determine the necessary angles to fully define the shape
% of the robot.  This script utilizes the equler lagrange approach

%CHANGE

clc
clear all
close all

%create symbolic variables
syms kappa2 kappa3 kappa4 L2 L3 L4 a1 t2 t3 t4 m2 m3 m4 g K3 Kdel g_dist t03 t0del

%% forward kineamtics
%determine bar lengths
n = 4;
K = [kappa2; kappa3; kappa4];
L = [L2; L3; L4];
t0 = [t2;t3;t4];
a(1)=a1;
for i = 1:n-1
    a(i+1) = (2*sin((L(i)*K(i))/2))/K(i);
end

x1_vec = [0;0;0];
x2_vec = [0;0;0];
x3_vec = [a(2)*cos(t2);
          a(2)*cos(t2);
          0];
x4_vec = [a(1);0;0];

%rotation matricies
R{1} = [1, 0, 0;
        0, 1, 0;
        0, 0, 1];
for i = 1:n-1
    R{i+1} = [cos(t0(i)), -sin(t0(i)), 0;
              sin(t0(i)),  cos(t0(i)), 0;
                     0,           0, 1];
end

% Determine COM location
C{1} = 0;
for i = 1:n-1
    phi = L(i)*K(i)*.5;
    A =(2*sin(phi))/(L(i)*K(i)^2);
    B = cos(phi)/K(i);
    C{i+1} = A-B;
end

%COM locations
P0G2 = a(2)/2*R{2}(1:3,1) - C{2}*R{2}(1:3,2);
P0G3 = a(3)/2*R{3}(1:3,1) - C{3}*R{3}(1:3,2);
P0G4 = a(4)/2*R{4}(1:3,1) + C{4}*R{4}(1:3,2);

%% Geometry
f_squared = a(1)^2+a(2)^2-2*a1*a(2)*cos(t2);
del = acos((a(3)^2+a(4)^2-f_squared)/(2*a(3)*a(4)));
g_dist = a(3)-a(4)*cos(del);
h = a(4)*sin(del);
r = a(1)-a(2)*cos(t2);
s = a(2)*sin(t2);
t3_final = atan((h*r-g_dist*s)/(g_dist*r+h*s))
t4_final = del+t3_final


%% deterine Y
Y = [                          1;
    (diff(t3_final,t2));
    (diff(t4_final,t2))];

%% Determine gravitational forces
%determine jacobains at COM
JVG2 = [(diff(P0G2,t2)), (diff(P0G2,t3)), (diff(P0G2,t4))];
JVG3 = [(diff(P0G3,t2)), (diff(P0G3,t3)), (diff(P0G3,t4))];
JVG4 = [(diff(P0G4,t2)), (diff(P0G4,t3)), (diff(P0G4,t4))];

%determine g0
g_dir = [0 -g 0].'; %gravity expressed as a 3-space vector
g0 = -[JVG2.' JVG3.' JVG4.']*[m2.*g_dir; m3.*g_dir; m4.*g_dir];

g_a = (subs(Y.'*g0,{t3,t4},{t3_final,t4_final}));

%% Deterime spring forces
tau0 = [            0;
    K3*(pi-t2+t3-t03);
    Kdel*(t4-t3-t0del)];

tau_a = (subs(Y.'*tau0,{t3,t4},{t3_final,t4_final}));

%% Find the residual
R = g_a-tau_a;
dR_dt2 = (diff(R,t2));








