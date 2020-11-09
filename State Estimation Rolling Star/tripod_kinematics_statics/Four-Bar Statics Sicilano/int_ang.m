function [thC] = int_ang(Pa,Pb,Pc)
%int_ang takes in 3 points and determines the internal angle between the
%lines AB and BC
%
%   Inputs: 
%
%   Pa = first point
%   Pb = second point
%   Pc = third point
%
%   Outputs:
%
%   thC = the angle beween the lines AB and BC, theta C
%   


% side lengths of triangle
a = norm(Pa-Pb);
b = norm(Pb-Pc);
c = norm(Pa-Pc);

%calculate angle
thC = acos((a^2+b^2-c^2)/(2*a*b));
end

