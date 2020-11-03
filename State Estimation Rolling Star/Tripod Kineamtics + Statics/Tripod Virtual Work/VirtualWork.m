% VirtualWork.m
% solves fourbar problem with springs at pins B and C

% Prepare Workspace
clear variables; close all; clc;

%% initialize linkage
a = 110;         % crank length (m)
b = 195;         % coupler length (m)
c = 120;         % rocker length (m)
d = 190;         % length between ground pins (m)


figure; hold on; grid on
%% initialize linkage and find starting configuration
test = true;
switch test
  case true
    axis equal
    q0 = 225*pi/180;
    [q,fval,~,output]  = fsolve(@(q)VWork(q,a,b,c,d),q0);  % use fsolve to find initial configuration
    [theta3,theta4,g,h] = fourbar(q,a,b,c,d);
    xA = [0;0];
    xB = a*[cos(q); sin(q)];
    xC = xB + b*[cos(theta3); sin(theta3)];
    xD = [d; 0];
    xD2 = xC - c*[cos(theta4); sin(theta4)];
    plot([xA(1) xB(1)],[xA(2) xB(2)],'LineWidth',2)
    plot([xB(1) xC(1)],[xB(2) xC(2)],'LineWidth',2)
    plot([xC(1) xD(1)],[xC(2) xD(2)],'LineWidth',2)
    plot([xC(1) xD2(1)],[xC(2) xD2(2)],'LineWidth',2)
  case false
    [theta2,Phi,J,dPhi] = deal(zeros(1,361));
    for i = 1:361
      theta2(i) = (i-1)*pi/180;
      [Phi(i),J(i)] = VWork(theta2(i),a,b,c,d);
    end
    dTheta2 = pi/180;
    for i = 2:361
      dPhi(i) = (Phi(i) - Phi(i-1))/dTheta2;
    end

 plot(theta2*180/pi,dPhi,'LineWidth',2)
 plot(theta2*180/pi,J,'.')
 % plot(theta2*180/pi,beta4,'LineWidth',1) 
 % plot(theta2*180/pi,dd4,'o','MarkerSize',8)
  %  plot(theta2*180/pi,delta3,'LineWidth',1) 
  %  plot(theta2*180/pi,delta4,'o','MarkerSize',8)
    xlim([0 360]);
end

function [W,J] = VWork(theta2,a,b,c,d)

%% mass and spring properties
  m2 = 0;  m3 = 0; m4 = 0.4; g = 9.81;
  F2 = [0; -m2*g]; F3 = [0; -m3*g]; F4 = [1; -m4*g];
  kB = 0; kC = 100; thetaB = -pi/2; thetaC = pi/2;
  
%% solve for theta3 and theta4 at current guess  
  [theta3,theta4,h,delta] = fourbar(theta2,a,b,c,d);
  phiB = theta3 - theta2 - thetaB;  phiC = theta4 - theta3 - thetaC;
  [e2,n2] = UnitVector(theta2); [e3,n3] = UnitVector(theta3); [e4,n4] = UnitVector(theta4);
  
%% partials w.r.t. theta2  
  delta3 = a*c*sin(theta2 - theta4)/(b*h);
  delta4 = a*sin(theta2 - theta3)/h;

  d31 = delta3 - 1;
  d43 = delta4 - delta3;
  
  C = d43/tan(delta);
  beta4 = delta4*((delta3-1)/tan(theta3-theta2) - C);
  beta3 = delta3*((delta4-1)/tan(theta4-theta2) - C);
  b43 = beta4 - beta3;
  
%% virtual work of forces
  Wf2 = dot(F2,a*n2/2);
  Wf3 = dot(F3,a*n2 + b*n3*delta3/2);
  Wf4 = dot(F4,c*n4*delta4/2);

%% virtual work of springs  
  WkB = -kB*phiB*d31;
  WkC = -kC*phiC*d43;

  W = Wf2 + Wf3 + Wf4 + WkB + WkC;
  
%% calculate Jacobian
  JkB = -kB*(d31^2 + phiB*beta3);
  JkC = -kC*(d43^2 + phiC*b43);
  Jf2 = dot(F2,a*e2/2);
  Jf3 = dot(F3,a*e2 + (b/2)*(e3*delta3^2 - n3*beta3));
  Jf4 = dot(F4,       (c/2)*(e4*delta4^2 - n4*beta4));
  J = JkB + JkC - Jf2 - Jf3 - Jf4;
  
end

function [e,n] = UnitVector(theta)
  e = [cos(theta); sin(theta)];
  n = [-e(2); e(1)];
end
function [theta3, theta4, h,delta] = fourbar(theta2,a,b,c,d)

  r = d - a*cos(theta2);
  s = a*sin(theta2);
  f2 = r^2 + s^2;                     % f squared
  delta = acos((b^2+c^2-f2)/(2*b*c)); % angle between coupler and rocker

  g = b - c*cos(delta);
  h = c*sin(delta);

  theta3 = atan2((h*r - g*s),(g*r + h*s));
  theta4 = theta3 + delta;

end