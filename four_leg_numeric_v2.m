%%RISS 2020 FWD KIN 4 LEGGED ROLLING STAR numeric%%

clc
clear all
close all 

%% enter constants for sample calculations
%enter arbitrary values for thetas
t1 = 0 * pi/180;
t1p = 60 * pi/180;

% determine bar length form curvature, k
L = .05;   %length of curved arm (m)
k1 = 10;   %curvature of arm 1 (m^-1)
k2 = 12;   %curvature of arm 2 (m^-1)
k3 = 10;   %curvature of arm 3 (m^-1)
k1p = 10;  %curvature of arm 1p (m^-1)
K = [k1,k2,k3,k1p]; %vector of curvatures
a = barcalc(K,L); %determine bar lengths (m)
%assign bar lengths to variables
a1 = a(1);
a2 = a(2);
a3 = a(3);
a1p = a(4);

%% closed loop constraints
x = sqrt(a1^2 + a1p^2 - 2*a1*a1p*cos(t1p-t1));
alpha = asin((a1p*sin(t1p-t1))/x);
beta = acos((x^2 + a2^2 - a3^2)/(2*x*a2));
t2 = pi - alpha - beta;

%% forward kinematics frame 0 to 2
a02 = a(1:2,1);
alpha = [0 0].';
d = [0 0].';
theta = [t1 t2].';
[T0_n,Tnm1_n] = fwdkinRISSnum(a02, alpha, d, theta);

%% forward kinematics fram 0 to 1p
a01p = a(4,1);
alpha = [0].';
d = [0].';
theta = [t1p].';
[T0_1p,T0_1] = fwdkinRISSnum(a01p, alpha, d, theta);

%% check other angles
t1 = t1*180/pi
t2 = t2*180/pi
t3 = (pi - asin((x*sin(beta))/a3))*180/pi
t1p = t1p * 180/pi

%% plot square
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

plot(x(1:2),y(1:2),'g-',x(2:3),y(2:3),'b-',x(3:4),y(3:4),'-',x(4:5),y(4:5),'r-')


%% detemine COM
for i = 1:length(x)-1
    m(i) = (y(i+1) - y(i))/(x(i+1) - x(i));       
    m_new = -1/m(i);
    xc = (x(i+1) + x(i))/2;
    yc = (y(i+1) + y(i))/2;
    alpha = (L*K(i))/2;    
    z = (sin(alpha))/(K(i)*alpha) - cos(alpha)/K(i);
    low_thresh = 1e-10;
    high_thresh = 1e10;
    if abs(m(i)) < low_thresh
        if i <= 2
            yf(i) = yc + z;
        else
            yf(i) = yc - z;
        end
        xf(i) = xc;
        
    elseif abs(m(i)) > high_thresh
        if i<=2 
            xf(i) = xc - z;
        else
            xf(i) = xc +z;
        end
        yf(i) = yc;
    else            
        if i<=2
            xf(i) = -sqrt((z^2)/(1+m_new^2))+xc;
        elseif i == 3 && m(i)<0
            xf(i) = -sqrt((z^2)/(1+m_new^2))+xc;
        else
            xf(i) = sqrt((z^2)/(1+m_new^2))+xc;
        end
        yf(i) = (xf(i)-xc)*m_new+yc;
    end
end
x_cent = sum(xf)/length(xf)
y_cent = sum(yf)/length(yf)

hold on
plot(xf,yf,'x',x_cent,y_cent,'ro')
legend('arm 1','arm 2','arm 3','arm 1p','COM of arm','COM of robot')
axis equal

%% check condition from book
t1 = t1*pi/180;
t2 = t2*pi/180;
t3 = t3*pi/180;
a03 = a(1:3,1);
alpha = [0 0 0].';
d = [0 0 0].';
theta = [t1 t2 t3].';
[T0_n,Tnm1_n] = fwdkinRISSnum(a03, alpha, d, theta);

P0_3 = T0_n{3}(1:3,4);
P0_1p = T0_1p{1}(1:3,4);
R3_0 = T0_n{3}(1:3,1:3).';

cond = vpa((R3_0*(P0_3-P0_1p)),4)



