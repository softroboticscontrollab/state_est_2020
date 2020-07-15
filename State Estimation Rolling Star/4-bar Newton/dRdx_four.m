function [J]=dRdx_four(z,param);
%
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

J(1,1)=-1;
J(1,2)=1;
J(1,3)=1;
J(1,4)=0;
J(1,5)=1;

J(2,1)=-a1;
J(2,2)=a1;
J(2,3)=0;
J(2,4)=-(a2*(K1*(a1^2+a1p^2-2*a1*a1p*cos(z(4)))^(3/2)*sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+K2*a1p*(-cos(z(4))*a1p+a1)*(a1*cos(z(4))-a1p))*sqrt(-(-2*a1*a1p*cos(z(4))+a1p^2+(a2+a3+a1)*(-a2-a3+a1))*(-2*a1*a1p*cos(z(4))+a1p^2+(a2-a3+a1)*(-a2+a3+a1))/((a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2^2))-sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))*K2*a1*a1p*sin(z(4))*(-2*a1*a1p*cos(z(4))+a1^2-a2^2+a3^2+a1p^2))/(sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))*(a1^2+a1p^2-2*a1*a1p*cos(z(4)))^(3/2)*sqrt(-(-2*a1*a1p*cos(z(4))+a1p^2+(a2+a3+a1)*(-a2-a3+a1))*(-2*a1*a1p*cos(z(4))+a1p^2+(a2-a3+a1)*(-a2+a3+a1))/((a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2^2))*a2);
J(2,5)=0;

J(3,1)=0;
J(3,2)=-1;
J(3,3)=0;
J(3,4)=-((a2*(cos(z(4))*a1p-a1)*(a1*cos(z(4))-a1p)*(m2+m3)*sin(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2)))-2*sin(z(4))*(a1*a1p*cos(z(4))-(1/2)*a1^2-(1/2)*a1p^2)*(m3+m1p)*sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))*sqrt(-(-2*a1*a1p*cos(z(4))+a1p^2+(a2+a3+a1)*(-a2-a3+a1))*(-2*a1*a1p*cos(z(4))+a1p^2+(a2-a3+a1)*(-a2+a3+a1))/((a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2^2))-2*sin(z(4))*(m2+m3)*(a1*a1p*cos(z(4))-(1/2)*a1^2+(1/2)*a2^2-(1/2)*a3^2-(1/2)*a1p^2)*sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))*a1*sin(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2))))*g*a1p/(2*sqrt(-(-2*a1*a1p*cos(z(4))+a1p^2+(a2+a3+a1)*(-a2-a3+a1))*(-2*a1*a1p*cos(z(4))+a1p^2+(a2-a3+a1)*(-a2+a3+a1))/((a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2^2))*(a1^2+a1p^2-2*a1*a1p*cos(z(4)))^(3/2)*sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))*a1);
J(3,5)=0;

J(4,1)=0;
J(4,2)=-1;
J(4,3)=-1;
J(4,4)=0;
J(4,5)=0;

J(5,1)=a1;
J(5,2)=0;
J(5,3)=0;
J(5,4)=-((g*a1p*a2*(-a1*cos(z(4))+a1p)*(cos(z(4))*a1p-a1)*(m2+m3)*sin(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2)))-(a1^2+a1p^2-2*a1*a1p*cos(z(4)))^(3/2)*(g*a1p*(m1+m3)*sin(z(4))+2*K1)*sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))-2*a1p*K2*(-a1*cos(z(4))+a1p)*(cos(z(4))*a1p-a1))*a2*sqrt(-(-2*a1*a1p*cos(z(4))+a1p^2+(a2+a3+a1)*(-a2-a3+a1))*(-2*a1*a1p*cos(z(4))+a1p^2+(a2-a3+a1)*(-a2+a3+a1))/((a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2^2))-(g*a2*(m2+m3)*sin(asin(a1p*sin(z(4))/sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4))))+acos((a2^2-a3^2+a1^2+a1p^2-2*a1*a1p*cos(z(4)))/(2*sqrt(a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2)))-2*K2)*a1p*sin(z(4))*a1*(-2*a1*a1p*cos(z(4))+a1^2-a2^2+a3^2+a1p^2)*sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4)))))/(2*sqrt((-cos(z(4))*a1p+a1)^2/(a1^2+a1p^2-2*a1*a1p*cos(z(4))))*sqrt(-(-2*a1*a1p*cos(z(4))+a1p^2+(a2+a3+a1)*(-a2-a3+a1))*(-2*a1*a1p*cos(z(4))+a1p^2+(a2-a3+a1)*(-a2+a3+a1))/((a1^2+a1p^2-2*a1*a1p*cos(z(4)))*a2^2))*(a1^2+a1p^2-2*a1*a1p*cos(z(4)))^(3/2)*a2);
J(5,5)=0;
end


