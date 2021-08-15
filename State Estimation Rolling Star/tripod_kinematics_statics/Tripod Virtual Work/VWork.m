function [W,J]=VWork(theta2,param)
%   Authors: Dr. Eric Constans (Rose-Hulman)
%            Sam Alvares (Rose-Hulman)
%   Date: 11/11/2020
%
%   Description: VWork solves for the virtual work and jacobain given
%   configuration parameters and the crank angle
%
%   Inputs:
%
%   theta2 = angle from ground to crank
%   param = vector of all parameters
%
%   Outputs:
%
%   W = residual, vidrual work
%   J = jacobian (1x1)

%% Assign parameters
%   param = [a b c d m2 m3 m4 g kB kC thetaB thetaC F(i) Fg2 Fg3 Fg4 gravity force_to_right];
a = param(1);
b = param(2);
c = param(3);
d = param(4);
m2 = param(5);
m3 = param(6);
m4 = param(7);
g = param(8);
kB = param(9);
kC = param(10);
thetaB = param(11);
thetaC = param(12);
F = param(13);
Fg2 = param(14);
Fg3 = param(15);
Fg4 = param(16);
gravity = param(17);
force_to_right = param(18);

%% Determine external forces
if (gravity == true) && (force_to_right == true)
    F2 = [F; Fg2]; F3 = [0; Fg3]; F4 = [0; Fg4];
elseif (gravity == true) && (force_to_right == false)
    F2 = [0; Fg2]; F3 = [0; Fg3]; F4 = [0; Fg4];
elseif (gravity == false) && (force_to_right == true)
    F2 = [F; 0]; F3 = [0; 0]; F4 = [0; 0];
else % gravity == false and force_to_right == false)
    F2 = [0; 0]; F3 = [0; 0]; F4 = [0; 0];
end

%% Solve for pose at current guess
[theta3,theta4,h,delta] = fourbar(theta2,a,b,c,d);
[e2,n2] = UnitVector(theta2); [e3,n3] = UnitVector(theta3); [e4,n4] = UnitVector(theta4);

%% Partials w.r.t. theta2
delta3 = a*c*sin(theta2 - theta4)/(b*h);
delta4 = a*sin(theta2 - theta3)/(h);

d31 = delta3 - 1;
d43 = delta4 - delta3;

C = d43/tan(delta);
beta4 = delta4*((delta3-1)/tan(theta3-theta2) - C);
if (theta4==theta2)
    beta3 = 0;
else
    beta3 = delta3*((delta4-1)/tan(theta4-theta2) - C);
end
b43 = beta4 - beta3;

%% Virtual work of forces
Wf2 = dot(F2,a*n2/2);
Wf3 = dot(F3,a*n2 + b*n3*delta3/2);
Wf4 = dot(F4,c*n4*delta4/2);

%% Virtual work of springs
% first determine the net deflection of the spring
phiB = theta3 - theta2 - thetaB;  phiC = theta4 - theta3 - thetaC;

%determine the virtual work due to the springs
WkB = -kB*phiB*d31;
WkC = -kC*phiC*d43;

%% Determine the net virtual work
W = Wf2 + Wf3 + Wf4 + WkB + WkC;

%% Calculate Jacobian
JkB = -kB*(d31^2 + phiB*beta3);
JkC = -kC*(d43^2 + phiC*b43);
Jf2 = dot(F2,a*e2/2);
Jf3 = dot(F3,a*e2 + (b/2)*(e3*delta3^2 - n3*beta3));
Jf4 = dot(F4,       (c/2)*(e4*delta4^2 - n4*beta4));
J = JkB + JkC - Jf2 - Jf3 - Jf4;

end