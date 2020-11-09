function Y = closure_Y_fourbar(param)
%closure_Y_fourbar Calculate the loop-closure Jacobian "Ypsilon" from the
%Siciliano text.
%
%   See Sam's other code for what needs to be in "params." It's a vector of
%   things to sub into the symbolically-solved expression.

%create variables for each of the parameters
L1=param(1);
L2=param(2);
L3=param(3);
a4=param(4);
kappa1 = param(5);
kappa2 = param(6);
kappa3 = param(7);
m1=param(8);
m2=param(9);
m3=param(10);
g=param(11);
K2=param(12);
K3=param(13);
t02=param(14);
t03=param(15);

% This is 3 x 1.

Y = [1;
    - (a4*(4*cos(t1) - 4*cos((L1*kappa1)/2)^2*cos(t1) + 2*a4*kappa1*sin((L1*kappa1)/2) + a4^2*kappa1^2*cos(t1) + 2*a4*kappa1*sin((L1*kappa1)/2)*cos(t1)^2))/(kappa1^2*(1 - (a4^2*kappa1^2*sin(t1)^2)/(a4^2*kappa1^2 + 4*cos(t1)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2))^(1/2)*((a4^2*kappa1^2 - 4*cos((L1*kappa1)/2)^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1) + 4)/kappa1^2)^(3/2)) - (a4*sin((L1*kappa1)/2)*sin(t1)*(4*kappa1^2*kappa2^2 - 4*kappa1^2*kappa3^2 + 4*kappa2^2*kappa3^2 - 4*kappa2^2*kappa3^2*cos((L1*kappa1)/2)^2 + 4*kappa1^2*kappa3^2*cos((L2*kappa2)/2)^2 - 4*kappa1^2*kappa2^2*cos((L3*kappa3)/2)^2 + a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2)*cos(t1)))/(2*kappa1^3*kappa2*kappa3^2*sin((L2*kappa2)/2)*((a4^2*kappa1^2 - 4*cos((L1*kappa1)/2)^2 + 4*a4*kappa1*sin((L1*kappa1)/2)*cos(t1) + 4)/kappa1^2)^(3/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2));
    -(a4*sin((L1*kappa1)/2)*sin(t1)*(4*kappa2^2*kappa3^2*sin((L1*kappa1)/2)^2 - 4*kappa1^2*kappa3^2*sin((L2*kappa2)/2)^2 - 4*kappa1^2*kappa2^2*sin((L3*kappa3)/2)^2 + a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2)*cos(t1)))/(kappa1^3*kappa3*sin((L2*kappa2)/2)^2*sin((L3*kappa3)/2)*((a4^2*kappa1^2 + 4*cos(t1)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)/kappa1^2)^(1/2)*(1 - (kappa1^2*kappa2^2*((4*sin((L1*kappa1)/2)^2)/kappa1^2 + (4*sin((L2*kappa2)/2)^2)/kappa2^2 - (4*sin((L3*kappa3)/2)^2)/kappa3^2 + a4^2 + (4*a4*sin((L1*kappa1)/2)*cos(t1))/kappa1)^2)/(16*sin((L2*kappa2)/2)^2*(a4^2*kappa1^2 + 4*cos(t1)*a4*kappa1*sin((L1*kappa1)/2) + 4*sin((L1*kappa1)/2)^2)))^(1/2)*((a4^2*kappa1^2*kappa2^2*kappa3^2 + 4*cos(t1)*a4*kappa1*kappa2^2*kappa3^2*sin((L1*kappa1)/2) - 4*kappa1^2*kappa2^2*sin((L3*kappa3)/2)^2 - 4*kappa1^2*kappa3^2*sin((L2*kappa2)/2)^2 + 4*kappa2^2*kappa3^2*sin((L1*kappa1)/2)^2)^2/(kappa1^4*kappa2^2*kappa3^2*sin((L2*kappa2)/2)^2*sin((L3*kappa3)/2)^2))^(1/2))];

end

