function [dR]=dresid_single_pend(x,param)
%   Inputs:
%
%   x = variable to solve for
%   param = vector of all parameters
% 
%   Outputs:
%
%   dR = slop of residual 

%create variables for each of the parameters
K=param(1);
m=param(2);
g=param(3);
l=param(4);
theta0=param(5);
%calculate the slope of the residual
dR=-K + m*g*l*sin(x);
end