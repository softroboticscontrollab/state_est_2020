%% RISS 2020 statics 4-bar numeric solver %%
% Using the solution of four_bar_statics_v5, find the pose of the 4-bar
% mechanism at static equilibrium

%% Prep the workspace
clc
clear all
close all

%% Determine setup configuration
COM_on_bar = true;
gravity = false;
neutral_right_angles = true;
force_to_right = false;
test_case = 2;
t1i=pi/1.8;  % initial guess of t1 (rad)

%% Enter parameters for numeric solver
n = 4;      % number bars
m1=1;       % mass limb 1(kg)
m2=1;       % mass limb 2(kg)
m3=1;       % mass limb 3(kg)
g=9.81;     % acceleration due to gravity (m/s/s)
k2=10;      % spring constant at joint 2 (m N/rad)
k3=10;      % spring constant at joint 3 (m N/rad)
F = 1;     % force to the right of the left bar. Only used if force_to_right == true
if (COM_on_bar)
    if (test_case == 1)
        a1 = 4;
        a2 = 4;
        a3 = 5;
        a4 = 2;
    elseif (test_case == 2)
        a1 = 4;
        a2 = 5;
        a3 = 3;
        a4 = 6;
    elseif (test_case == 3)
        % Simplified energy verification
        a1 = 1;
        a2 = 1;
        a3 = 1;
        a4 = 1;
    end
    % Create a vector of bar lengths
    a = [a1, a2, a3];
else
    L1 = 5;     % arc leangth limb 1 (m)
    L2 = 5;     % arc leangth limb 2 (m)
    L3 = 5;     % arc leangth limb 3 (m)
    a4 = 1;     % ground bar length (m)
    kappa1 = .9;    % curvature limb 1 (1/m)
    kappa2 = .4;    % curvature limb 2 (1/m)
    kappa3 = .1;    % curvature limb 3 (1/m)
    % determine bar lengths
    K = [kappa1; kappa2; kappa3]; % create vector of curvature
    L = [L1; L2; L3];             % create vector of arc lengths
    a = zeros(n-1,1);
    for i = 1:n-1
        a(i) = barcalc(K(i),L(i)); % caluclate virtual bar length
    end
end

%calculate initial spring lengths at t2i, the initial guess for t2
% x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1i));
% alpha = asin((a4*sin(t1i))/x);
% beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
% zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
% t2i = pi - alpha - beta;  % initial t2 (rad)
% t3i = pi - zeta;          % initial t3 (rad)

if (neutral_right_angles)
    t02 = pi/2;          % unstretched spring length for t2 (rad)
    t03 = pi/2;
else
    t02 = 1;          % unstretched spring length for t2 from calibration (rad)
    t03 = 1;          % unstretched spring length for t3  from calibration (rad)
end

% Create a vector of parameters for the numeric solver
if (COM_on_bar)
    param = [a1 a2 a3 a4 m1 m2 m3 g k2 k3 t02 t03 F force_to_right gravity COM_on_bar]; %create a vector of all paramters
else
    param = [L1 L2 L3 a4 kappa1 kappa2 kappa3 m1 m2 m3 g k2 k3 t02 t03 gravity COM_on_bar]; %create a vector of all paramters
end

%% Run Newton Raphson numeric solver
R=@resid_four_bar_statics;      % residual function
dRdx=@dresid_four_bar_statics;  % slope of residual
tol=0.0000001;                  % solution tolerance
maxIter=1000000;                % max iteratios to find solution
toggle=1;                       % 1 prints guesses

% Run numeric solver
[t1_sol,er_est,t1_guess]=func_newton(R,dRdx,t1i,tol,maxIter,toggle,param);

%% Forward Kinematics
% find the numeric value of t2 and t3 using the output of the numeric
% solver
x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1_sol));
alpha = asin((a4*sin(t1_sol))/x);
beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
t2_sol = pi - alpha - beta; % settling angle for joint 2 (rad)
t3_sol = pi - zeta;         % settling angle for joint 3 (rad)

%print joint angles
fprintf('\ntheta1: %1.4f\ntheta2: %1.4f\ntheta3: %1.4f\n',t1_sol,t2_sol,t3_sol)

% Determine forward kinamtics at the final pose
theta = [t1_sol;t2_sol; t3_sol];
zero = [0;0;0];
[T0_n,Tnm1_n] = fwdkinRISSnum(a, zero, zero, theta);

% Create corner vector at final position
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

% Determine forward kinamtics at the inital pose
% theta_initial = [t1i;t2i;t3i];
% zero = [0;0;0];
% [T0_n_initial,Tnm1_n_initial] = fwdkinRISSnum(a, zero, zero, theta_initial);

% Create corner vector for initial pose
% xin(1) = 0;
% yin(1) = 0;
% xin(2) = a(1)*cos(t1i);
% yin(2) = a(1)*sin(t1i);
% xin(3) = xin(2) + a(2)*cos(t1i+t2i);
% yin(3) = yin(2) + a(2)*sin(t1i+t2i);
% xin(4) = xin(3) + a(3)*cos(t1i+t2i+t3i);
% yin(4) = yin(3) + a(3)*sin(t1i++t2i+t3i);
% xin(5) = 0;
% yin(5) = 0;

%% Find COM Locations
% Determine the location of the COM at the final pose
P0Gx = zeros(n-1,1);
P0Gy = zeros(n-1,1);
for i = 1:n-1
    if (COM_on_bar)
        C = 0;
    else
        phi = L(i)*K(i)*.5;
        A =(2*sin(phi))/(L(i)*K(i)^2);
        B = cos(phi)/K(i);
        C = A-B;
    end
    if i==1
        P0Gx(1) = a(1)/2*T0_n{1}(1,1)+C*T0_n{1}(1,2);
        P0Gy(1) = a(1)/2*T0_n{1}(2,1)+C*T0_n{1}(2,2);
    else
        P0Gx(i) = T0_n{i-1}(1,4)+a(i)/2*T0_n{i}(1,1)+C*T0_n{i}(1,2);
        P0Gy(i) = T0_n{i-1}(2,4)+a(i)/2*T0_n{i}(2,1)+C*T0_n{i}(2,2);
    end
end

% Determine the location of the COM at the initial pose
% P0Gx_initial = zeros(n-1,1);
% P0Gy_initial = zeros(n-1,1);
% for i = 1:n-1
%     if (COM_on_bar)
%         C = 0;
%     else
%         phi = L(i)*K(i)*.5;
%         A =(2*sin(phi))/(L(i)*K(i)^2);
%         B = cos(phi)/K(i);
%         C = A-B;
%     end
%     if i==1
%         P0Gx_initial(1) = a(1)/2*T0_n_initial{1}(1,1)+C*T0_n_initial{1}(1,2);
%         P0Gy_initial(1) = a(1)/2*T0_n_initial{1}(2,1)+C*T0_n_initial{1}(2,2);
%     else
%         P0Gx_initial(i) = T0_n_initial{i-1}(1,4)+a(i)/2*T0_n_initial{i}(1,1)+C*T0_n_initial{i}(1,2);
%         P0Gy_initial(i) = T0_n_initial{i-1}(2,4)+a(i)/2*T0_n_initial{i}(2,1)+C*T0_n_initial{i}(2,2);
%     end
% end

%% Simple energy verification
if (force_to_right == true) && (test_case == 3) && (gravity == false)
    % Code is aplicable when force_to_right == true and test_case == 3 and gravity == false
    % Check energy for simplified case
    % phi2 = pi/2-t2_sol;
    % phi3 = pi/2-t3_sol;
    % Ei = g*((m1+m2)*.5*a1+m3*a1);
    % Ef = g*((m1+m2)*.5*a1+m3*a1)*sin(t1_sol)+.5*10*(phi2^2+phi3^2);
    % fprintf('\nInitiail energy: %.9f\n', Ei)
    % fprintf('Final energy: %.9f\n', Ef)
    
    % Print results for simplified energy verification (case 3)
    g0 = [ 10*a2*sin(t1_sol + t2_sol) + 10*a1*sin(t1_sol) + 5*a3*sin(t1_sol + t2_sol + t3_sol);
        10*a2*sin(t1_sol + t2_sol) + 5*a3*sin(t1_sol + t2_sol + t3_sol);
        5*a3*sin(t1_sol + t2_sol + t3_sol)]
    
    ga = 10*a2*sin(asin((a4*sin(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) - t1_sol + acos((a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)))) + 10*a1*sin(t1_sol) - (10*a2*sin(asin((a4*sin(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) - t1_sol + acos((a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)))) - 5*a3*sin(t1_sol - asin((a4*sin(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) + acos((a1^2 + 2*cos(t1_sol)*a1*a4 - a2^2 - a3^2 + a4^2)/(2*a2*a3)) - acos((a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)))))*(((a1*a4*sin(t1_sol))/(a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) - (a1*a4*sin(t1_sol)*(a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2))/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(3/2)))/(1 - (a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)))^(1/2) + ((a4*cos(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2) + (a1*a4^2*sin(t1_sol)^2)/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(3/2))/(1 - (a4^2*sin(t1_sol)^2)/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2))^(1/2)) - 5*a3*sin(t1_sol - asin((a4*sin(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) + acos((a1^2 + 2*cos(t1_sol)*a1*a4 - a2^2 - a3^2 + a4^2)/(2*a2*a3)) - acos((a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)))) - (5*a1*a4*sin(t1_sol - asin((a4*sin(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) + acos((a1^2 + 2*cos(t1_sol)*a1*a4 - a2^2 - a3^2 + a4^2)/(2*a2*a3)) - acos((a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2))))*sin(t1_sol))/(a2*(1 - (a1^2 + 2*cos(t1_sol)*a1*a4 - a2^2 - a3^2 + a4^2)^2/(4*a2^2*a3^2))^(1/2))
    
    
    tau0=[                   0;
        -k2*(t02 + t2_sol - pi);
        -k3*(t03 + t3_sol - pi)]
    
    taua = - k2*(((a1*a4*sin(t1_sol))/(a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) - (a1*a4*sin(t1_sol)*(a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2))/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(3/2)))/(1 - (a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)))^(1/2) + ((a4*cos(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2) + (a1*a4^2*sin(t1_sol)^2)/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(3/2))/(1 - (a4^2*sin(t1_sol)^2)/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2))^(1/2))*(asin((a4*sin(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) - t02 + acos((a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)))) - (a1*a4*k3*sin(t1_sol)*(t03 - pi + acos((a1^2 + 2*cos(t1_sol)*a1*a4 - a2^2 - a3^2 + a4^2)/(2*a2*a3))))/(a2*a3*(1 - (a1^2 + 2*cos(t1_sol)*a1*a4 - a2^2 - a3^2 + a4^2)^2/(4*a2^2*a3^2))^(1/2))
    
    
    
    Y = [
        
    1;
    - ((a1*a4*sin(t1_sol))/(a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2)) - (a1*a4*sin(t1_sol)*(a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2))/(2*a2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(3/2)))/(1 - (a1^2 + 2*cos(t1_sol)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)))^(1/2) - ((a4*cos(t1_sol))/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(1/2) + (a1*a4^2*sin(t1_sol)^2)/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2)^(3/2))/(1 - (a4^2*sin(t1_sol)^2)/(a1^2 + 2*cos(t1_sol)*a1*a4 + a4^2))^(1/2);
    (a1*a4*sin(t1_sol))/(a2*a3*(1 - (a1^2 + 2*cos(t1_sol)*a1*a4 - a2^2 - a3^2 + a4^2)^2/(4*a2^2*a3^2))^(1/2))]

% word and spring energy for simple energy verification
phi2 = pi/2-t2_sol;
phi3 = pi/2-t3_sol;
W = 10*a(3)/2*cos(t3_sol)
Ef = .5*10*(phi2^2+phi3^2)

end

%% Plot
% Plot COM at initial guess
% figure(1)
% initial_COM = plot(P0Gx_initial,P0Gy_initial,'kx');

% Plot pose at initial guess
% figure(1)
% initial_corners = plot(xin,yin,'--k');

% Plot COM at final pose
hold on
final_COM = plot(P0Gx,P0Gy,'rx');

% Plot final pose
hold on
final_corners = plot(x,y,'r');

% Create legend
% legend('Inital COM','Final COM','Initial Pose','Final Pose')
legend('Final COM','Final Pose')
xlabel('x-position (m)')
ylabel('y-position (m)')
axis equal

%% Verify solution with energy
% Determine initial PE with unstretched spring lengths
% t1_guess = t1_guess(1:length(t1_guess)-1);
% for i = 1:length(t1_guess)
%     %find the numeric value of t2 and t3
%     x = sqrt(a(1)^2 + a4^2 - 2*a(1)*a4*cos(pi-t1_guess(i)));
%     alpha = asin((a4*sin(t1_guess(i)))/x);
%     beta = acos((x^2 + a(2)^2 - a(3)^2)/(2*x*a(2)));
%     % zeta = asin((x*sin(beta))/a(3));
%     zeta = acos((-(x)^2+a(2)^2+a(3)^2)/(2*a(2)*a(3)));
%     t2_guess = pi - alpha - beta;
%     t3_guess = pi - zeta;
%     % t3_sol = atan2(s_zeta,c_zeta);
%
%     %determine forward kinamtics
%     theta = [t1_guess(i);t2_guess; t3_guess];
%     zero = [0;0;0];
%     [T0_n,Tnm1_n] = fwdkinRISSnum(a, zero, zero, theta);
%
%     for j = 1:n-1
%         phi = L(j)*K(j)*.5;
%         A =(2*sin(phi))/(L(j)*K(j)^2);
%         B = cos(phi)/K(j);
%         C = A-B;
%         if j==1
%             P0Gx(1) = a(1)/2*T0_n{1}(1,1)+C*T0_n{1}(1,2);
%             P0Gy(1) = a(1)/2*T0_n{1}(2,1)+C*T0_n{1}(2,2);
%         else
%             P0Gx(j) = T0_n{j-1}(1,4)+a(j)/2*T0_n{j}(1,1)+C*T0_n{j}(1,2);
%             P0Gy(j) = T0_n{j-1}(2,4)+a(j)/2*T0_n{j}(2,1)+C*T0_n{j}(2,2);
%         end
%     end
%
%     Ugrav(i) = g*(m1*P0Gy(1)+m2*P0Gy(2)+m3*P0Gy(3));
%     Uspring(i) = .5*k2*(pi-t2_guess-t02)^2+.5*(pi-t3_guess-t03)^2;
%     Utot(i) = Ugrav(i) + Uspring(i);
% end
% Ugrav
% Uspring
% Utot
%
% x_axis = 1:length(t1_guess);
% figure
% subplot(3,1,1);
% plot(x_axis,Ugrav,'b-')
% ylabel('Ugrav')
% hold on
% subplot(3,1,2)
% plot(x_axis,Uspring,'k-')
% ylabel('Uspring')
% hold on
% subplot(3,1,3)
% plot(x_axis,Utot,'r-')
% ylabel('Utot')
%
t2_constans = 2*pi-(pi-t1_sol)-(pi-t2_sol)-(pi-t3_sol)





