function [soln,er_est,x]=func_newton(R,dRdx,xi,tol,maxIter,toggle,param)
%   Inputs:
%
%   R = function to deterrmine the residuals
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
%set x as initial guess
x(i)=xi;
%numeric solver
while i<maxIter
    I(i)=i;
    Rxi(i)=R(x(i),param);
    dRdxi(i)=dRdx(x(i),param);
    corr(i)=abs((Rxi(i)/dRdxi(i)));
    x(i+1)=x(i)-(Rxi(i)/dRdxi(i)); 
    if corr(i) < tol
        soln=x(i);
        er_est=corr(i);
        i=maxIter+1;
    else   
    i=i+1;
    end
end
%print solution if toggle is on
 if toggle == 1 
        fprintf(' Count       xi        R(xi)      dRdx(xi)       corr      \n');
        for y=1:length(I)
        fprintf('  %3.0f %10.3f %12.3f %11.3f %12.3f  \n',I(y), x(y), Rxi(y), dRdxi(y), corr(y));
        
        end
 end
    
    