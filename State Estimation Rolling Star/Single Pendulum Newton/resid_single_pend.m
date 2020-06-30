function [R]=resid_single_pend(x,param)
K=param(1);
m=param(2);
g=param(3);
l=param(4);
theta0=param(5);
R=K*(theta0-x)-m*g*l*cos(x);
end