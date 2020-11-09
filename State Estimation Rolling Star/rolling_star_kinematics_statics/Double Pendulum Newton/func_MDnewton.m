function [xrt,er_est]=func_MDnewton(resid_vec,dRdx,xi,tol,maxIter,toggle,param)
%   Inputs:
%
%   resid_vec = function to deterrmine vector of residuals
%   dRdx = rate of change of residuals
%   xi = solver initial value
%   tol = solution tolerance
%   maxIter = max interations to find solution
%   toggle = 1 prints results
%   param = vector of parameters
% 
%   Outputs:
%
%   soln = solution
%   er_est = error

%initialize i
i=1;
%set solution to initial guess
x=xi;
%run numeric solver
if toggle == 1
fprintf('\nIteration      corr       xrt\n');
end
if rank(dRdx(x,param)) < length(xi)
    fprintf('choose new initial guess');
    xrt=0;
    er_est=0;
else
while i<maxIter
    Rxi=resid_vec(x,param)';
    dRdxi=dRdx(x,param);
    corr=inv(dRdxi)*(-Rxi);
     if toggle == 1
        g=sprintf('   % .2e ',x);
        fprintf('  %3f  % 11.2e%s\n',i, max(abs(corr)),g);
     end
       x=x+corr;
    if max(abs(corr)) < tol
        xrt=x;
        er_est=corr;
        i=maxIter+1;
    else   
    i=i+1;
    end
end
end
end


    
    