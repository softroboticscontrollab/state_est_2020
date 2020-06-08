%%RISS 2020 STATICS EXAMPLE%%

clc
clear all
close all 

%create variables
syms a1 a2 t1 t2 m1 m2 g k1 k2 L C

%% forward kinematics fram 0 to 2
a = [a1 a2].';
alpha = [0 0].';
d = [0 0].';
theta = [t1 t2].';
[T0_n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta);

%% Determine Jacobians at COM
%enter necessary directions and vectors
x0_1 = T0_n{1}(1:3,1);
y0_1 = T0_n{1}(1:3,2);
x0_2 = T0_n{2}(1:3,1);
y0_2 = T0_n{2}(1:3,2);
P0_1 = T0_n{1}(1:3,4);
z0_0 = [0 0 1].';
z0_1 = [0 0 1].';

%calculate zeta, linar distance from chord to COM
K = [k1 k2];
for i = 1:length(K)
    alpha = (L*K(i))/2;    
    zeta(i) = (sin(alpha))/(K(i)*alpha) - cos(alpha)/K(i);
end

P0_G1 = (a1/2)*x0_1+zeta(1)*y0_1;
P0_G2 = P0_1 + (a2/2)*x0_2+zeta(2)*y0_2;
P1_G2 = (a2/2)*x0_2+zeta(2)*y0_2;

%for link 1
JV1 = simplify(cross(z0_0,P0_G1));
JV2 = [0 0 0].';

%create J matrix
JVG1 = [JV1 JV2];

%For link 2
JV1 = simplify(cross(z0_0,P0_G2));
JV2 = simplify(cross(z0_1,P1_G2));

JVG2 = [JV1 JV2];

%% create gravity vector
g_dir = [0 -g 0].';
G = -[JVG1.' JVG2.']*[m1.*g_dir; m2.*g_dir];

%% determine restoring force
F = C*[t1 - 51.43*pi/180; t2 - 51.43*pi/180]

%% Determine numeric solution
% determine bar length form curvature, k
L_num = .05;   %length of curved arm (m)
k1_num = 20;   %curvature of arm 1 (m^-1)
k2_num = 60;   %curvature of arm 2 (m^-1)
K = [k1_num,k2_num]; %create a vector of curvatures (m^-1)
a = barcalc(K,L_num); %determine bar lengths (m)

%enter other constants
C_num = 10; %torsional spring constant N-m/rad
a1_num = a(1); %bar 1 length (m^-1)
a2_num = a(2); %bar 2 length (m^-1)
m1_num = 10;   %mass bar 1 (kg)
m2_num = 20;   %mass bar 2 (kg)
g_num = 9.81;  %gravity (kg-m/s/s)
t1_num = 30 * pi/180; %theta 1 (rad)
t2_num = 60 * pi/180; %theta 2 (rad)

%determine the numeric gravity vector
G_num = vpa(subs(G,{L,k1,k2,a1,a2,m1,m2,g,t1,t2},{L_num,k1_num,k2_num,a1_num,a2_num,m1_num,m2_num,g_num,t1_num,t2_num}),4)

%determine restoring force numerically
F_num = vpa(subs(F,{C,t1,t2},{C_num,t1_num,t2_num}),6)

%calculate joint torques
T_num = vpa(G_num+F_num,6)

%detereine where joints are and plot
T0_1_num = vpa(subs(T0_n{1},{a1,a2,t1,t2},{a1_num,a2_num,t1_num,t2_num}),4);
T0_2_num = vpa(subs(T0_n{2},{a1,a2,t1,t2},{a1_num,a2_num,t1_num,t2_num}),4);
x(1) = 0;
y(1) = 0;
x(2) = T0_1_num(1,4);
y(2) = T0_1_num(2,4);
x(3) = T0_2_num(1,4);
y(3) = T0_2_num(2,4);

plot(x(1:2),y(1:2),'g-',x(2:3),y(2:3),'b-')

%detemine COM of links and plot
P0_G1_num = vpa(subs(P0_G1,{a1,t1,k1,L},{a1_num,t1_num,k1_num,L_num}),6);
P0_G2_num = vpa(subs(P0_G2,{a1,a2,t1,t2,k2,L},{a1_num,a2_num,t1_num,t2_num,k2_num,L_num}),6);

xCOM = [0 P0_G1_num(1,1) P0_G2_num(1,1)];
yCOM = [0 P0_G1_num(2,1) P0_G2_num(2,1)];

hold on
plot(xCOM,yCOM,'x')
legend('arm 1','arm 2','COM of arm')

