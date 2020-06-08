% This function takes in vectors of k, curvature, and L, arc length and 
% calculates the bar length
function a = barcalc(K,L)
% intialize T
a = zeros(length(K),1);
for i=1:length(a)
   a(i) = (2*sin((L*K(i))/2))/K(i);
end

end