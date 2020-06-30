function [J]=dRdx(x,param);
K1 = param(1);
K2 = param(2);
m1 = param(3);
m2 = param(4);
g = param(5);
L1 = param(6);
L2 = param(7);
theta0=param(8);

J(1,1)=K1-m1*g*(L1/2)*sin(x(1))-m2*g*L1*sin(x(1))-m2*g*(L2/2)*sin(x(1)+x(2));
J(1,2)=-m2*g*(L2/2)*sin(x(2));
J(2,1)=0;
J(2,2)=K2-m2*g*(L2/2)*sin(x(1)+x(2));
end


