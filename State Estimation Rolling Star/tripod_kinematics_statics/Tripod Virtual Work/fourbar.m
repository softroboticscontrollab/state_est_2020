function [theta3, theta4, h, delta] = fourbar(theta2,a,b,c,d)
%   Authors: Dr. Eric Constans (Rose-Hulman)
%            Sam Alvares (Rose-Hulman)
%   Date: 11/11/2020
%
%   Description: fourbar detemirmines the the geometry of a fourbar
%   mechanisim given all of the bar lengths and the joint angle of the
%   crank. This program uses the computationally efficient aproach found in
%   Introduction to Mechanisms - With Computer Applications 
%
%   Inputs:
%
%   theta2 = the angle from the ground to the crank
%   a-d = bar lengths of the mechinism
%
%   Outputs:
%
%   thetai = joint angle
%   h = geometry output
%   

  r = d - a*cos(theta2);
  s = a*sin(theta2);
  f2 = r^2 + s^2;                     % f squared
  delta = acos((b^2+c^2-f2)/(2*b*c)); % angle between coupler and rocker

  g = b - c*cos(delta);
  h = c*sin(delta);
  
  theta3 = atan2((h*r - g*s),(g*r + h*s));
  theta4 = theta3 + delta;

end