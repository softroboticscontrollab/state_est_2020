function [e,n] = UnitVector(theta)
%   Authors: Dr. Eric Constans (Rose-Hulman)
%            Sam Alvares (Rose-Hulman)
%   Date: 11/11/2020
%
%   Description: determine the unit vector given a joint angle
%
%   Inputs:
%
%   theta = joint angle
%
%   Outputs:
%
%   e = unit vector 
%   n = unit normal vector
%   
  e = [cos(theta); sin(theta)];
  n = [-e(2); e(1)];
end