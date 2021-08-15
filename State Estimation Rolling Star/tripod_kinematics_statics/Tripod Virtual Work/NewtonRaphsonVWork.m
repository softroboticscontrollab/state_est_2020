function [soln]=NewtonRaphsonVWork(VWork,theta2_0,tol,maxIter,toggle,param)
%   Authors: Dr. Eric Constans (Rose-Hulman)
%            Sam Alvares (Rose-Hulman)
%   Date: 11/11/2020
%
%   Description: NewtonRaphsonVWork is a Newton Raphson solver that allows
%   the pose of the robot to be determined given an initial guess fof the 
%   angle.
%
%   Inputs:
%
%   VWork = function to determine the virtual work and jacobian at a
%           particular pose
%   theta2_0 = solver initial value (rad)
%   tol = solution tolerance
%   maxIter = max interations to find solution
%   toggle = 1 prints results
%   param = vector of parameters
%
%   Outputs:
%
%   soln = solution for theta2 (rad)
%   er_est = error

%% Initialize
i=1;
%set theta2 as initial guess
theta2(i)=theta2_0;

%% Run numeric solver
while i<maxIter
    I(i)=i; %vector of count number
    [Wi(i),Ji(i)]=VWork(theta2(i),param); % determine W and J at current guess
    corr(i)=abs((Wi(i)/Ji(i))); % correction for next guess
    theta2(i+1)=theta2(i)-(Wi(i)/Ji(i));
    if corr(i) < tol
        soln=theta2(i);
        er_est=corr(i);
        i=maxIter+1;
    else
        i=i+1;
    end
end

%% Print solution if toggle is on
if toggle == 1
    fprintf(' Count    theta2i    W(theta2i)  J(theta2i)      corr      \n');
    for y=1:length(I)
        fprintf('  %3.0f %10.3f %12.3f %11.3f %12.3f  \n',I(y), theta2(y), Wi(y), Ji(y), corr(y));
    end
end

