%% RISS 2020 4-Bar Robot Statics %%
% Given the arc lengths, curvatures, material proprties and necessary 
% angle to fully define the robot's pose, determine the reaction forces 
% at the frames of the base of the robot

clc
clear all
close all 

%% Enter knowns (arbitrary) 
% curvatures (m^-1)
k1 = 15;
k2 = 15;
k3 = 15;
k1p = 15;
% arc length (m)
L = .05;
% angles
t1p = 110 * (pi/180);
t1 = 0;
% acceleration due to gravity (m/s/s)
g=9.81;
% masses (kg)
m1 = .3;
m2 = .3;
m3 = .3;
m1p = .3;

K = [k1,k2,k3,k1p]; %vector of curvatures
a = barcalc(K,L); %determine bar lengths (m)

%assign bar lengths to variables
a1 = a(1);
a2 = a(2);
a3 = a(3);
a1p = a(4);

%% Closed loop constraints
x = sqrt(a1^2 + a1p^2 - 2*a1*a1p*cos(t1p-t1));
alpha = asin((a1p*sin(t1p-t1))/x);
beta = acos((x^2 + a2^2 - a3^2)/(2*x*a2));
%calculate unknown angle (note: t3 not needed in future calcs)
t2 = pi - alpha - beta;

%% Forward kinematics frame 0 to 2
a02 = a(1:2,1);
alpha = [0 0].';
d = [0 0].';
theta = [t1 t2].';
[T0_n,Tnm1_n] = fwdkinRISSnum(a02, alpha, d, theta);

x0 = 0;
y0 = 0;
x1 = T0_n{1}(1,4);
y1 = T0_n{1}(2,4);
x2 = T0_n{2}(1,4);
y2 = T0_n{2}(2,4);
%% Forward kinematics fram 0 to 1p
[T0_1p,T0_1] = fwdkinRISSnum(a(4,1), 0, 0, t1p);

x1p = T0_1p{1}(1,4);
y1p = T0_1p{1}(1,4);

%% Statics
f2y = g*((1/2)*m1*a1+(1/2)*m1p*cos(t1p)*a1p+m3*((x1p+x2)*(1/2))+m2*((x1+x2)*(1/2)))/a1;
f1y = -f2y+g*(m1+m2+m3+m1p);

%% Plot bars 
x(1) = 0;
y(1) = 0;
for i = 1:length(T0_n)
    x(i+1) = T0_n{i}(1,4);
    y(i+1) = T0_n{i}(2,4);
end
x(i+2) = T0_1p{1}(1,4);
y(i+2) = T0_1p{1}(2,4);
x(i+3) = 0;
y(i+3) = 0;

plot(x(1:2),y(1:2),'g-',x(2:3),y(2:3),'k',x(3:4),y(3:4),'m',x(4:5),y(4:5),'c')

%% Detemine and plot COM's
for i = 1:length(x)-1
    xf(i) = (x(i)+x(i+1))/2;
    yf(i) = (y(i)+y(i+1))/2;
end
x_cent = sum(xf)/length(xf);
y_cent = sum(yf)/length(yf);

hold on
plot(xf,yf,'x',x_cent,y_cent,'ro')


%% Plot force
scale = 200; %scale force so it will fit on plot
hold on
quiver(0,0,0,f1y/scale,'b-')
hold on
quiver(x1,0,0,f2y/scale,'r-')

legend('arm 1','arm 2','arm 3','arm 1p','COM of arm','COM of robot','f1y','f2y')
axis equal
xlabel('x-location (m)')
ylabel('y-location (m)')
set(gcf,'color','w');
axis ([-.04 .07 0 .07])

