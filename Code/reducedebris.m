function [ deltaX ] = reducedebris(X, reductionConstant, N)
%reducedebris 
%   Reduces the debris in each layer according so ome constant

numlevels = size(X,1)/2;

numLasersPerLevel = zeros(numlevels,1);
counter = 0;
currIndex = numlevels;
cutoffIndex = 2;
if numlevels <= 2
   cutoffIndex = 1;
end
    
while counter < N
   if currIndex == cutoffIndex
      currIndex = numlevels;
   end
   numLasersPerLevel(currIndex) = numLasersPerLevel(currIndex) + 1;
   
   currIndex = currIndex - 1;
   counter = counter + 1;
end

deltaX(2:2:numlevels*2) = numLasersPerLevel .* reductionConstant;

end

