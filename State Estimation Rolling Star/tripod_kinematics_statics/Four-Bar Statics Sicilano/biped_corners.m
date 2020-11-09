function corners = biped_corners(vertices_limbs, vertices_feet)
%biped_corners Calculate the joint positions for the biped (four-bar)
%   
%   Given vertex data (vertices), this script outputs a 2 x 4 matrix of the
%   (x,y) coordinates of each of the four joints of the virtual four-bar
%   mechanism.
%
%   Inputs:
%       vertices_limbs, indexed as vertices_limbs{j}{x, :) etc
%       vertices_feet, indexed as vertices_feet{j}(x, :) etc
%

% NOTE HERE that the limb ordering has been swapped with respect to the DER
% simulation: limb 1 is the rightmost, limb 2 is top, limb 3 is left.
% Foot 1 is right foot, foot 2 is left foot.

% We assume the bottom joints to be at the geometric mean in the
% x-direction for the feet, and that the other two corners are at the first
% vertices of limbs 2 and 1. See Sam's diagram, it helps.

% Node 1 is bottom right foot:
% First foot, nodes 1 and 2, x coordinate
x1 = (1/2) * (vertices_feet{1}(1,2) - vertices_feet{1}(1,1));
x1 = x1 + vertices_feet{1}(1,1);
% second foot, same
x4 = (1/2) * (vertices_feet{2}(1,2) - vertices_feet{2}(1,1));
x4 = x4 + vertices_feet{2}(1,1);

% top right corner = first limb node 1
x2 = vertices_limbs{1}(1,1);
y2 = vertices_limbs{1}(2,1);
% top left corner = second limb node 1
x3 = vertices_limbs{2}(1,1);
y3 = vertices_limbs{2}(2,1);


corners = [ x1, 0;
            x2, y2;
            x3, y3;
            x4, 0]';
end

