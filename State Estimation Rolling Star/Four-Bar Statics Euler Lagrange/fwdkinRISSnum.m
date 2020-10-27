function [T0n,Tnm1_n] = fwdkinRISSnum(a, alpha, d, theta)
%   Inputs:
%
%   a, alpha, d, theta = vectors of numeric DH parameters
% 
%   Outputs:
%   
%   T0n = structure of all numeric tranfromation matricies from frame 0 to
%         frame n
%   Tnm1_n = structure of all numeric tranfromation matricies from frame n
%            minus 1 to the next consecutive fram, frame n
%

% intialize T structures as a structure of 4x4 identity matricies
for j = 1:length(a)
    T{j} = eye(4);  %intermedite transformations from fram 0 to n
    T_int{j} = eye(4); %intermedite transformations from fram n-1 to n
end

% determine the first entry in both T and T_int
T{1} = [cos(theta(1))  -sin(theta(1))*cos(alpha(1))     sin(theta(1))*sin(alpha(1))     a(1)*cos(theta(1));
          sin(theta(1))   cos(theta(1))*cos(alpha(1))    -cos(theta(1))*sin(alpha(1))     a(1)*sin(theta(1));
          0               sin(alpha(1))                   cos(alpha(1))                   d(1);
          0               0                               0                               1];

T_int{1} = [cos(theta(1))  -sin(theta(1))*cos(alpha(1))     sin(theta(1))*sin(alpha(1))     a(1)*cos(theta(1));
          sin(theta(1))   cos(theta(1))*cos(alpha(1))    -cos(theta(1))*sin(alpha(1))     a(1)*sin(theta(1));
          0               sin(alpha(1))                   cos(alpha(1))                   d(1);
          0               0                               0                               1];
      
% determine the rest of the transformations        
for i=2:length(a)
    T_int{i} = [cos(theta(i))  -sin(theta(i))*cos(alpha(i))     sin(theta(i))*sin(alpha(i))     a(i)*cos(theta(i));
              sin(theta(i))   cos(theta(i))*cos(alpha(i))    -cos(theta(i))*sin(alpha(i))     a(i)*sin(theta(i));
              0               sin(alpha(i))                   cos(alpha(i))                   d(i);
              0               0                               0                               1];
    T{i}=(T{i-1}*T_int{i});
end

% return T0n and Tnm1_n
T0n = T;
Tnm1_n = T_int;
end