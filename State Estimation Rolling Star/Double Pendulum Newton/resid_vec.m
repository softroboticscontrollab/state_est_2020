function [Rvec]=resid_vec(x,param)
%   Inputs:
%
%   x = vector of variables to solve for
%   param = vector of all parameters
% 
%   Outputs:
%
%   Rvec = vector of slop of residuals

%create variables for each of the parameters
K1 = param(1);
K2 = param(2);
m1 = param(3);
m2 = param(4);
g = param(5);
L1 = param(6);
L2 = param(7);
theta0=param(8);
%calculate residuals 
Rvec(1)=K1*(x(1)-theta0)+m1*g*(L1/2)*cos(x(1))+m2*g*(L2/2)*cos(x(1)+x(2));
Rvec(2)=K2*(x(2)-theta0)+m2*g*(L2/2)*cos(x(1)+x(2));
end