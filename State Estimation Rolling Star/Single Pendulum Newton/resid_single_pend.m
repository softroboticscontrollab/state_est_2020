function [R]=resid_single_pend(x,param)
%   Inputs:
%
%   x = variable to solve for
%   param = vector of all parameters
% 
%   Outputs:
%
%   R = residual

%create variables for each of the parameters
K=param(1);
m=param(2);
g=param(3);
l=param(4);
theta0=param(5);

%caluclate residual
R=K*(theta0-x)-m*g*l*cos(x);
end