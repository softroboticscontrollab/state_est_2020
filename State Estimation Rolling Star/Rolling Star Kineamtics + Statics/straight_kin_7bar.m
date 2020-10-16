function [straight_tips,ms,Pc] = straight_kin_7bar(T0_n,T0_np,a,K,Ls)
%
%   straight_kin_7bar calculates the location of the tips of the stright
%   sections of the rolling star robot knowing the forward kineamtics and
%   geometry
%
%   Inputs:
%
%   T0_n = Transformation matricies from 0 to joint n
%
%   T0_np = Transformation matrix from 0 to joint 1'
%
%   a = stright bar lengths [a1,a2,a3,a4,a5,a6,a1']
%
%   K = arm curvatures [k1,k2,k3,k4,k5,k6,k1']
%
%   Ls = length of stright section
%
%   Outputs:
%
%   stright_tips = 2x7 matrix of xy coordinates of the tips of the robot.
%   R1 = x, R2 = y
%
%   ms = bisecting slope of the two curves
%
%   Pc = center of curvatures for all arms, struct of 3x1 vectors
%

% define constants
num_limbs = 7;

%% determine curvature centers of all arms
% intialize Pc structures as a 3x1 vectors of zeros
for i = 1:num_limbs
    Pc{i} = zeros(3,1);
end

% calculate curvature centers.  Separate case for arm 1 bc base jt is
% assumed to 0,0.
for i = 1:num_limbs-1
    if i == 1
        Pc{i} = (a(i)/2)*T0_n{i}(1:3,1)-sqrt(-(a(i)/2)^2+(1/K(i))^2)*T0_n{i}(1:3,2);
    else
        Pc{i} = T0_n{i-1}(1:3,4)+(a(i)/2)*T0_n{i}(1:3,1)-sqrt(-(a(i)/2)^2+(1/K(i))^2)*T0_n{i}(1:3,2);
    end
end
% calculate curvatu8re center for arm 1'
Pc{num_limbs} = (a(num_limbs)/2)*T0_np{1}(1:3,1)+sqrt(-(a(num_limbs)/2)^2+(1/K(num_limbs))^2)*T0_np{1}(1:3,2);

%% Determine the slopes of all the straight sections
% Initialize the bisecting slope as a vector of zeros
ms = zeros(num_limbs,1);
% Calculate bisecting curvatures
for i = 1:num_limbs
    if i < num_limbs
        % slope of the n arm at joint n
        mn = (Pc{i}(1)-T0_n{i}(1,4))/(T0_n{i}(2,4)-Pc{i}(2));
        % slope of the n+1 arm at joint n
        mnp1 = (Pc{i+1}(1)-T0_n{i}(1,4))/(T0_n{i}(2,4)-Pc{i+1}(2));
        % calculate bisecting slope - if slopes are opposite signs, use
        % correct case
        if (mn < 0 && mnp1 > 0)
            ms(i) = (mn*mnp1-1+sqrt((mn^2+1)*(mnp1^2+1)))/(mn+mnp1);
        elseif (mn > 0 && mnp1 < 0)
            ms(i) = (mn*mnp1-1-sqrt((mn^2+1)*(mnp1^2+1)))/(mn+mnp1);
        else
            %https://demonstrations.wolfram.com/AngleBisectorsOfTwoIntersectingLines/
            ms(i) = (mn*mnp1-1+sqrt((mn^2+1)*(mnp1^2+1)))/(mn+mnp1);
        end
        
    else %for 1' arm
        % slope of the first arm at base jt
        m1 = Pc{1}(1)/-Pc{1}(2);
        % slope of the 1' arm at base jt
        m1p = Pc{i}(1)/-Pc{i}(2);
        if (m1 < 0 && m1p > 0) || (m1 > 0 && m1p < 0)
            ms(i) = (abs(m1)+abs(m1p))/2;
        else
            ms(i) = (m1*m1p-1+sqrt((m1^2+1)*(m1p^2+1)))/(m1+m1p);
        end
    end
end

%% Calcuate tip locations
% initialize matrix for xy straight tip locations as zeros
straight_tips = zeros(2,7);
% Calcuate tip locations
for i = 1:num_limbs-1
    A = sqrt((Ls)^2/(1+ms(i)^2));
    % tip location option 1
    P0t1_0 = [T0_n{i}(1,4) + A;
        T0_n{i}(2,4) + ms(i)*A;
        0];
    % tip location option 2
    P0t2_0 = [T0_n{i}(1,4) - A;
        T0_n{i}(2,4) - ms(i)*A;
        0];
    % calcuate the vector from jt n to the tip location in fram 0
    Pnt1_0 = P0t1_0 - T0_n{i}(1:3,4);
    % rotation matrix from from 0 to n
    R0_n = T0_n{i}(1:3,1:3).';
    % vector from jt n to the tip location in fram n
    Pnt1_1 = R0_n*Pnt1_0;
    % need tip location outside of enclosed shape.  If the y-component of
    % Pnt1_1 is negative, P0t1_0 is the correct solution, else the other
    % option,P0t2_0, is the solution
    if Pnt1_1(2)<0
        straight_tips(1,i) = P0t1_0(1);
        straight_tips(2,i) = P0t1_0(2);
    elseif Pnt1_1(2)>0
        straight_tips(1,i) = P0t2_0(1);
        straight_tips(2,i) = P0t2_0(2);
    else %Pnt1_1(2) == 0
        straight_tips(1,i) = P0t1_0(1);
        straight_tips(2,i) = P0t1_0(2);
    end
end

% for final stright section, similar code as above, but no need to
% transfrom to the nth frame, sincve the final stright section protrudes
% from the base fram
A = sqrt((Ls)^2/(1+(ms(num_limbs))^2));
P0t1_0 = [A;
    ms(num_limbs)*A;
    0];

P0t2_0 = [-A;
    -ms(num_limbs)*A;
    0];

if P0t1_0(2) < 0
    straight_tips(1,num_limbs) = P0t1_0(1);
    straight_tips(2,num_limbs) = P0t1_0(2);
elseif P0t1_0(2) > 0
    straight_tips(1,num_limbs) = P0t2_0(1);
    straight_tips(2,num_limbs) = P0t2_0(2);
else %P0t1_0(2) == 0
    straight_tips(1,num_limbs) = P0t1_0(1);
    straight_tips(2,num_limbs) = P0t1_0(2);
end

end