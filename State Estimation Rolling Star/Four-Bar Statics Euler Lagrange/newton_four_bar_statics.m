%% RISS 2020 statics 4-bar numeric solver %%
% Using the solution of four_bar_statics_v4, find the pose of the 4-bar
% mechanism at static equilibrium

%prep the workspace
clc
clear all
close all

%% enter parameters
n = 4; %number bars
L1 = 5;
L2 = 5;
L3 = 5;
a4 = 5;
kappa1 = .9;
kappa2 = .4;
kappa3 = .01;
m1=1;
m2=1;
m3=1;
g=9.81;
k2=10;
k3=10;

% determine bar lengths
K = [kappa1; kappa2; kappa3];
L = [L1; L2; L3];
for i = 1:n-1
    a(i) = barcalc(K(i),L(i));
end

%calculate initial spring lengths at t2 initial guess
t1i=pi/2.1;
x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1i));
alpha = asin((a4*sin(t1i))/x);
beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
% zeta = asin((x*sin(beta))/a(3));
zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
t02 = pi -( pi - alpha - beta);
t03 = pi - (pi - zeta);

param = [L1 L2 L3 a4 kappa1 kappa2 kappa3 m1 m2 m3 g k2 k3 t02 t03]; %create a vector of all paramters

%% newton raphson
%inputs and functions
R=@resid_four_bar_statics; %residual function
dRdx=@dresid_four_bar_statics; %slope of residual
% t1i=pi/2.1; %initial guess for solution
tol=0.00001; %solution tolerance
maxIter=1000; %max iteratios to find solution
toggle=1; %1 prints guesses

%Run numeric solver
[t1_sol,er_est,t1_guess]=func_newton(R,dRdx,t1i,tol,maxIter,toggle,param);

%% plot results
%find the numeric value of t2 and t3
x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1_sol));
alpha = asin((a4*sin(t1_sol))/x);
beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
% zeta = asin((x*sin(beta))/a(3));
zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
t2_sol = pi - alpha - beta;
t3_sol = pi - zeta;
% t3_sol = atan2(s_zeta,c_zeta);

%determine forward kinamtics
theta = [t1_sol;t2_sol; t3_sol];
zero = [0;0;0];
[T0_n,Tnm1_n] = fwdkinRISSnum(a, zero, zero, theta);

% Create corner vector
x(1) = 0;
y(1) = 0;
x(2) = a(1)*cos(t1_sol);
y(2) = a(1)*sin(t1_sol);
x(3) = x(2) + a(2)*cos(t1_sol+t2_sol);
y(3) = y(2) + a(2)*sin(t1_sol+t2_sol);
x(4) = x(3) + a(3)*cos(t1_sol+t2_sol+t3_sol);
y(4) = y(3) + a(3)*sin(t1_sol+t2_sol+t3_sol);
x(5) = 0;
y(5) = 0;

plot(x,y)
axis equal
axis([-7 2 -2 7])

%% Verify solution with energy
%determine the location of the COM
for i = 1:n-1
    phi = L(i)*K(i)*.5;
    A =(2*sin(phi))/(L(i)*K(i)^2);
    B = cos(phi)/K(i);
    C = A-B;
    if i==1
        P0Gx(1) = a(1)/2*T0_n{1}(1,1)+C*T0_n{1}(1,2);
        P0Gy(1) = a(1)/2*T0_n{1}(2,1)+C*T0_n{1}(2,2);
    else
        P0Gx(i) = T0_n{i-1}(1,4)+a(i)/2*T0_n{i}(1,1)+C*T0_n{i}(1,2);
        P0Gy(i) = T0_n{i-1}(2,4)+a(i)/2*T0_n{i}(2,1)+C*T0_n{i}(2,2);
    end
end
hold on
plot(P0Gx,P0Gy,'rx')


% determine initial pose
xi(1) = 0;
yi(1) = 0;
xi(2) = L(1)*cos(pi/2);
yi(2) = L(1)*sin(pi/2);
xi(3) = xi(2) + L(2)*cos(pi);
yi(3) = yi(2) + L(2)*sin(pi);
xi(4) = -a4;
yi(4) = 0;
xi(5) = 0;
yi(5) = 0;

% determine inital COM pos
P0Gxi = [0;-L2/2;-L2];
P0Gyi = [L1/2;L1;L1/2];

% Plot inital pose
hold on
plot(xi,yi,'--k')

%Plot initial COM position
hold on
plot(P0Gxi,P0Gyi,'kx')



% Determine initial PE with unstretched spring lengths
t1_guess = t1_guess(1:length(t1_guess)-1);
for i = 1:length(t1_guess)
    %find the numeric value of t2 and t3
    x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1_guess(i)));
    alpha = asin((a4*sin(t1_guess(i)))/x);
    beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
    % zeta = asin((x*sin(beta))/a(3));
    zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
    t2_guess = pi - alpha - beta;
    t3_guess = pi - zeta;
    % t3_sol = atan2(s_zeta,c_zeta);
    
    %determine forward kinamtics
    theta = [t1_guess(i);t2_guess; t3_guess];
    zero = [0;0;0];
    [T0_n,Tnm1_n] = fwdkinRISSnum(a, zero, zero, theta);
    
    for j = 1:n-1
        phi = L(j)*K(j)*.5;
        A =(2*sin(phi))/(L(j)*K(j)^2);
        B = cos(phi)/K(j);
        C = A-B;
        if j==1
            P0Gx(1) = a(1)/2*T0_n{1}(1,1)+C*T0_n{1}(1,2);
            P0Gy(1) = a(1)/2*T0_n{1}(2,1)+C*T0_n{1}(2,2);
        else
            P0Gx(j) = T0_n{j-1}(1,4)+a(j)/2*T0_n{j}(1,1)+C*T0_n{j}(1,2);
            P0Gy(j) = T0_n{j-1}(2,4)+a(j)/2*T0_n{j}(2,1)+C*T0_n{j}(2,2);
        end
    end
    
    Ugrav(i) = g*(m1*P0Gy(1)+m2*P0Gy(2)+m3*P0Gy(3));
    Uspring(i) = .5*k2*(pi-t2_guess-t02)^2+.5*(pi-t3_guess-t03)^2;
    Utot(i) = Ugrav(i) + Uspring(i);
end
Ugrav
Uspring
Utot

x_axis = 1:length(t1_guess);
figure
subplot(3,1,1);
plot(x_axis,Ugrav,'b-')
ylabel('Ugrav')
hold on
subplot(3,1,2)
plot(x_axis,Uspring,'k-')
ylabel('Uspring')
hold on
subplot(3,1,3)
plot(x_axis,Utot,'r-')
ylabel('Utot')





