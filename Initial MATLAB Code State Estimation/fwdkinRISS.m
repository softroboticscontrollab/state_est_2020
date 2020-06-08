% This function takes in vectors of a, alpha, d and theta and computes
% the forward kinematics for a robot arm with symbolic results
function [T0n,Tnm1_n] = fwdkinRISS(a, alpha, d, theta)
% intialize T
for j = 1:length(a)
    T{j} = eye(4);
    T_int{j} = eye(4);
end

T{1} = [cos(theta(1))  -sin(theta(1))*cos(alpha(1))     sin(theta(1))*sin(alpha(1))     a(1)*cos(theta(1));
          sin(theta(1))   cos(theta(1))*cos(alpha(1))    -cos(theta(1))*sin(alpha(1))     a(1)*sin(theta(1));
          0               sin(alpha(1))                   cos(alpha(1))                   d(1);
          0               0                               0                               1];

T_int{1} = [cos(theta(1))  -sin(theta(1))*cos(alpha(1))     sin(theta(1))*sin(alpha(1))     a(1)*cos(theta(1));
          sin(theta(1))   cos(theta(1))*cos(alpha(1))    -cos(theta(1))*sin(alpha(1))     a(1)*sin(theta(1));
          0               sin(alpha(1))                   cos(alpha(1))                   d(1);
          0               0                               0                               1];
      
      
for i=2:length(a)
    %determine transformation matrix
    Tnm1_n = [cos(theta(i))  -sin(theta(i))*cos(alpha(i))     sin(theta(i))*sin(alpha(i))     a(i)*cos(theta(i));
              sin(theta(i))   cos(theta(i))*cos(alpha(i))    -cos(theta(i))*sin(alpha(i))     a(i)*sin(theta(i));
              0               sin(alpha(i))                   cos(alpha(i))                   d(i);
              0               0                               0                               1];
    T_int{i} = Tnm1_n;
    T{i}=simplify(T{i-1}*Tnm1_n);
end

% return T0n
T0n = T;
Tnm1_n = T_int;
end