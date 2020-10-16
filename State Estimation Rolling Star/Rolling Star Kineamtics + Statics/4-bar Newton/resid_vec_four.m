function [Rvec]=resid_vec_four(z,param)
K1 = param(1);
K2 = param(2);
m1 = param(3);
m2 = param(4);
m3 = param(5);
m1p = param(6);
g = param(7);
a1 = param(8);
a2 = param(9);
a3 = param(10);
a1p = param(11);
theta0=param(12);

Rvec(1)=-g*m1-z(1)+z(2)+z(3)+z(5);
Rvec(2)=-K2*(t0-Pi+asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2)))+K1*(t0-z(4))-(1/2)*m1*g*a1+(z(2)-z(1))*a1;
Rvec(3)=-z(2)+g*((1/2)*m1*a1+(1/2)*m1p*a1p*cos(z(4))+m3*((1/2)*cos(z(4))*a1p+(1/2)*a1-(1/2)*a2*cos(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2))))+m2*(a1-(1/2)*a2*cos(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2)))))/a1;
Rvec(4)=-z(3)+g*(m1+m1p+m2+m3)-z(2);
Rvec(5)=-(1/2)*m1*g*cos(z(4))*a1p-m3*g*((1/2)*cos(z(4))*a1p+(1/2)*a1-(1/2)*a2*cos(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2))))-m2*g*(a1-(1/2)*a2*cos(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2))))+z(1)*a1-K1*(t0-z(4))+K2*(t0-Pi+asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2)));

end