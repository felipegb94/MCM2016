function [ history, currIndex, currDecayIndex, deltaX, change ] = decaylevels(history, X, currIndex, currDecayIndex, currTime, dt)
%decaylevels Summary of this function goes here
%   Detailed explanation goes here

numlevels = size(history,1);
deltaX = zeros(2*numlevels,1);
deltaXIndeces = 2:2:2*numlevels;
numlevels = size(history,1);
change = 0;
if ((currTime + dt) >= (currDecayIndex(numlevels)))
    numIndecesToAdvance = floor((currTime + dt)) - (currDecayIndex(numlevels)-1);
    for j = 1:numIndecesToAdvance
        for i = numlevels:-1:1
            deltaXIndex = deltaXIndeces(i);
            numDecayObjects = history(i,currDecayIndex(i));
            deltaX(deltaXIndex) = deltaX(deltaXIndex) - numDecayObjects; 
            if(i > 1)
                history(i-1,currIndex) = history(i-1,currIndex) + numDecayObjects;
                deltaX(deltaXIndex - 2) = deltaX(deltaXIndex - 2) + numDecayObjects;
            end

            
        end
        change = 1;
        currIndex = currIndex + 1;
        currDecayIndex = currDecayIndex + 1;
    end
end

end

