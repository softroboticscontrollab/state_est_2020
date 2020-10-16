function [R]=resid_four_bar_statics(t1,param)
%   Inputs:
%
%   x = variable to solve for
%   param = vector of all parameters
% 
%   Outputs:
%
%   R = residual

%create variables for each of the parameters
a1=param(1);
a2=param(2);
a3=param(3);
a4=param(4);
m1=param(5);
m2=param(6);
m3=param(7);
g=param(8);
k2=param(9);
k3=param(10);
t0=param(11);

%caluclate residual
R=((a3*g*m3*cos(asin((a4*sin(t1))/(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)) - t1 + acos((a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2))) + asin(((1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2))/a3)))/2 + (a2*g*m2*cos(asin((a4*sin(t1))/(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)) - t1 + acos((a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)))))/2)*((a4*(a1^2*cos(t1) + a1*a4*cos(t1)^2 + a1*a4 + a4^2*cos(t1)))/(((a1 + a4*cos(t1))^2/(a1^2 + 2*cos(t1)*a1*a4 + a4^2))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(3/2)) + (a1*a4*sin(t1)*(a1^2 + 2*cos(t1)*a1*a4 - a2^2 + a3^2 + a4^2))/(2*a2*(1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(3/2))) - g*m2*((a2*cos(asin((a4*sin(t1))/(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)) - t1 + acos((a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)))))/2 - a1*cos(t1)) + k2*((a4*(a1^2*cos(t1) + a1*a4*cos(t1)^2 + a1*a4 + a4^2*cos(t1)))/(((a1 + a4*cos(t1))^2/(a1^2 + 2*cos(t1)*a1*a4 + a4^2))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(3/2)) + (a1*a4*sin(t1)*(a1^2 + 2*cos(t1)*a1*a4 - a2^2 + a3^2 + a4^2))/(2*a2*(1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(3/2)))*(asin((a4*sin(t1))/(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)) - t0 + acos((a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)))) + (a1*g*m1*cos(t1))/2 - (a3*g*m3*cos(asin((a4*sin(t1))/(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)) - t1 + acos((a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2))) + asin(((1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2))/a3)))/2 + (a1*a4*g*m3*cos(asin((a4*sin(t1))/(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)) - t1 + acos((a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)/(2*a2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2))) + asin(((1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2))/a3))*sin(t1)*(a1^2 + 2*cos(t1)*a1*a4 - a2^2 - a3^2 + a4^2))/(2*a2^2*(1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*((a1^2 + 2*cos(t1)*a1*a4 - a2^2 - a3^2 + a4^2)^2/(a2^2*a3^2))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2)) - (a1*a4*k3*sin(t1)*(t0 - asin(((1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2))/a3))*(a1^2 + 2*cos(t1)*a1*a4 - a2^2 - a3^2 + a4^2))/(a2^2*a3*(1 - (a1^2 + 2*cos(t1)*a1*a4 + a2^2 - a3^2 + a4^2)^2/(4*a2^2*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)))^(1/2)*((a1^2 + 2*cos(t1)*a1*a4 - a2^2 - a3^2 + a4^2)^2/(a2^2*a3^2))^(1/2)*(a1^2 + 2*cos(t1)*a1*a4 + a4^2)^(1/2));
end