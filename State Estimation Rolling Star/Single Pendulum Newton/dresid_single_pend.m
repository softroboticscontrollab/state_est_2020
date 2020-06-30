function [dR]=dresid_single_pend(x,param)
K=param(1);
m=param(2);
g=param(3);
l=param(4);
theta0=param(5);
dR=-K + m*g*l*sin(x);
end