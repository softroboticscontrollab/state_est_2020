function [xrt,er_est]=func_MDnewton(resid_vec,dRdx,xi,tol,maxIter,toggle,param)
i=1;
x=xi;
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


    
    