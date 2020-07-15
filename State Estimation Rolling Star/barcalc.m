function a = barcalc(K,L)
%   Inputs:
%
%   K = vectors of k, curvature
%   L = arc length of curved section of a single bar
% 
%   Outputs:
%   
%   a = vector of bar lengths assuming constant curvature
%
a = zeros(length(K),1);
for i=1:length(a)
   a(i) = (2*sin((L*K(i))/2))/K(i);
end

end