function [soln,er_est]=func_newton(R,dRdx,xi,tol,maxIter,toggle,param)
i=1;
x(i)=xi;
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
 if toggle == 1 
        fprintf(' Count       xi        R(xi)      dRdx(xi)       corr      \n');
        for y=1:length(I)
        fprintf('  %3.0f %10.3f %12.3f %11.3f %12.3f  \n',I(y), x(y), Rxi(y), dRdxi(y), corr(y));
        
        end
 end
    
    