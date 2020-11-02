%% RISS 2020 4-bar method of IC %%
% 

% prep the workspace
clc
clear all
close all

%% Find differential kinematics symboliclly - Long Method
fprintf('**********LONG METHOD**********')

%create symbolic variables
syms a1 a2 a3 a4 t2_sym 

% x1_vec = [0;0;0];
% x2_vec = [0;0;0];
% x3_vec = [a(2)*cos(t2);
%           a(2)*cos(t2);
%           0];
% x4_vec = [a(1);0;0];

% Geometry
f_squared = a1^2+a2^2-2*a1*a2*cos(t2_sym);
del = acos((a3^2+a4^2-f_squared)/(2*a3*a4));
g_dist = a3-a4*cos(del);
h = a4*sin(del);
r = a1-a2*cos(t2_sym);
s = a2*sin(t2_sym);
t3_sym = atan((h*r-g_dist*s)/(g_dist*r+h*s));
t4_sym = del+t3_sym;

% Find differential kinematics
dt3_dt2_sym = diff(t3_sym,t2_sym)
dt4_dt2_sym = diff(t4_sym,t2_sym)

%% Solution from symbolic approach
% enter paramteres
a1_num = 1;
a2_num = 3;
a3_num = 3.5;
a4_num = 5;
a_num = [a1,a2,a3,a4];
t2 = pi/3;

% determine numeric angles
f_squared = a1_num^2+a2_num^2-2*a1_num*a2_num*cos(t2);
del = acos((a3_num^2+a4_num^2-f_squared)/(2*a3_num*a4_num));
g_dist = a3_num-a4_num*cos(del);
h = a4_num*sin(del);
r = a1_num-a2_num*cos(t2);
s = a2_num*sin(t2);
t3 = atan((h*r-g_dist*s)/(g_dist*r+h*s));
t4 = del+t3;

params = [a1_num,a2_num,a3_num,a4_num,t2,t3,t4];

% Determine the solution of the symbolic differential kineamtics
dt3_dt2_sym_sol = vpa(subs(dt3_dt2_sym,{a1,a2,a3,a4,t2_sym,t3_sym,t4_sym},{params}),4)
dt4_dt2_sym_sol = vpa(subs(dt4_dt2_sym,{a1,a2,a3,a4,t2_sym,t3_sym,t4_sym},{params}),4)

%% Find differential kinematics symboliclly - IC Method
fprintf('**********IC METHOD**********')
A = [0;0];
B = [a2*cos(t2_sym);
     a2*sin(t2_sym)];
C = [a1+a4*cos(t4_sym);
       a4*sin(t4_sym)];
D = [a1;
      0];

m2 = B(2)/B(1);
m3 = (C(2)-B(2))/(C(1)-B(1));
m4 = C(2)/(C(1)-a1);

I24 = [(B(1) - B(2)/m3);
        0];

I13 =[(a1*m4)/(m4-m2);
        (m2*a1*m4)/(m4-m2)];

alpha = -I24(1);

beta = norm(I13-B);

dt3_dt2_IC = -a2/beta
dt4_dt2_IC = alpha/(alpha+a1)

dt3_dt2_IC_sol = vpa(subs(dt3_dt2_IC,{a1,a2,a3,a4,t2_sym,t3_sym,t4_sym},{params}),4)
dt4_dt2_IC_sol = vpa(subs(dt4_dt2_IC,{a1,a2,a3,a4,t2_sym,t3_sym,t4_sym},{params}),4)

I24_sol = vpa(subs(I24,{a1,a2,a3,a4,t2_sym,t3_sym,t4_sym},{params}),4);
I13_sol = vpa(subs(I13,{a1,a2,a3,a4,t2_sym,t3_sym,t4_sym},{params}),4);


%% Plot result
% plot pose
x(1) = 0;
y(1) = 0;
x(2) = a2_num*cos(t2);
y(2) = a2_num*sin(t2);
x(3) = a1_num+a4_num*cos(t4);
y(3) = a4_num*sin(t4);
x(4) = a1_num;
y(4) = 0;
x(5) = 0;
y(5) = 0;

plot(x,y)
hold on
plot(I24_sol(1),I24_sol(2),'rx')
hold on 
plot(I13_sol(1),I13_sol(2),'bx')
axis equal

% plot IC



