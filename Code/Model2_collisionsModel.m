function [] = Model2_collisionsModel(numlevels,tend)

close all;
% numlevels = 3;
% tend = 6000;
currTime = 0;
[history, decaytimes, offset, altitudes] = initlevels(numlevels,tend);

currDecayIndex = offset + 1 - decaytimes;
currIndex = offset + 1;

counter = 1;
population = zeros(2*numlevels,5*tend);
time = zeros(1,5*tend);
alldt = rand(1,5*tend);

% prepare stochastic variables
X = zeros(2*numlevels,tend);
X(1:2:end,1) = GetInitialSatellites(numlevels);
X(2:2:end,1) = sum(history,2);
% create reaction vectors
rxns = zeros(2*numlevels,numlevels+1);
% fragment-fragment reactions
fragdebris = 50;
fraginds = sub2ind(size(rxns),2:2:2*numlevels,1:numlevels);
rxns(fraginds) = fragdebris;
%other reactions

% initialize randoms
rs = rand(1,1000*tend);
Pk = log(1./rand(size(rxns,2),1));
Tk = zeros(size(rxns,2),1);
currTimeBin = 1;
while currTime < tend-1
    disp(currTime)
    % call stochastic reactions function
    rnd = rs(counter);
    [change,delt,Pk,Tk] =...
        stochevent(X(:,currTimeBin),rxns,Pk,Tk,rnd);
    
    history(:,currIndex) = history(:,currIndex)+change(2:2:end);
    %calculate decay
    if delt>=1
        delt=.999;
    end
    [history, currIndex, currDecayIndex,decay,flag] = decaylevels_sub(history, currIndex, currDecayIndex, currTime, delt);
    
    currTime = currTime + delt;
    currTimeBin = currDecayIndex(numlevels);

    if flag
        X(:,currTimeBin) = X(:,currTimeBin-1) + decay;
    end
    X(:,currTimeBin) = X(:,currTimeBin) + change;
    
    population(:,counter) = X(:,currTimeBin);
    time(counter) = currTime;
    counter = counter + 1;
end

figure(1);
for plotind = 1:numlevels
subplot(numlevels,1,plotind);
semilogy(time(1:counter-1), population(2*plotind,1:counter-1), 'LineWidth',1)
t = sprintf('Level %d: %d - %d Km',plotind,altitudes(plotind),altitudes(plotind+1));
title(t);
xlabel('Time (days)');
ylabel('Total Debris');


% subplot(numlevels,2,2*plotind);
% plot(time(1:counter-1), population(2*plotind-1,1:counter-1), 'LineWidth',1)
% t = sprintf('Level %d: %d - %d Km',plotind,altitudes(plotind),altitudes(plotind+1));
% title(t);
% xlabel('Time (days)');
% ylabel('Total Spacecraft');
% inds = find(collt(2,:)==plotind);
% hold on
% plot(collt(1,inds),zeros(1,length(inds)),'r*')
% hold off
end
% function return



end

%%
function [ history, currIndex, currDecayIndex, deltaX, change ] = decaylevels_sub(history, currIndex, currDecayIndex, currTime, dt)
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
