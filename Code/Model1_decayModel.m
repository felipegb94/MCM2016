function [] = Model1_decayModel(numlevels,tend)
close all;
%numlevels = 3;
%tend = 6000;
currTime = 0;
[history, decaytimes, offset, altitudes] = initlevels(numlevels,tend);

currDecayIndex = offset + 1 - decaytimes;
currIndex = offset + 1;

counter = 1;
population = zeros(numlevels,5*tend);
time = zeros(1,5*tend);
alldt = rand(1,5*tend);
X = zeros(2*numlevels,tend);
initPopulation = sum(history,2);
X(2:2:end,1) = initPopulation;

while currTime < tend
    if(mod(floor(currTime),1000) == 0)
        disp(currTime)
    end
    dt = alldt(counter);
    [history, currIndex, currDecayIndex, deltaX, change] = decaylevels(history, X, currIndex, currDecayIndex, currTime, dt);
    currTime = currTime + dt;
    currTimeBin = currDecayIndex(numlevels);
    if(change)
       X(:,currTimeBin) = X(:,currTimeBin-1) + deltaX;
    end
    population(:,counter) = X(2:2:end,currTimeBin);
    time(counter) = currTime;
    counter = counter + 1;
end

% Plotting %
%%%%%%%%%%%%
figure(1);
for plotind = 1:numlevels
subplot(numlevels,1,plotind);
semilogy(time(1:counter-1), population(plotind,1:counter-1), 'LineWidth',1)
t = sprintf('Level %d: %d - %d Km',plotind,altitudes(plotind),altitudes(plotind+1));
title(t);
xlabel('Time (days)');
ylabel('Total Debris');
xlim([0 tend])


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

% figure;
% subplot(3,1,1);
% semilogy(time(1:counter-1), population(1,1:counter-1), 'LineWidth',1)
% t = strcat('First Level (',num2str(altitudes(1)),'-',num2str(altitudes(2)),'km)');
% title(t);
% xlabel('Time (days)');
% ylabel('Total Debris');
% subplot(3,1,2);
% semilogy(time(1:counter-1), population(2,1:counter-1),'LineWidth',1)
% t = strcat('Second Level (',num2str(altitudes(2)),'-',num2str(altitudes(3)),'km)');
% title(t);
% xlabel('Time (days)');
% ylabel('Total Debris');
% subplot(3,1,3);
% plot(time(1:counter-1), population(3,1:counter-1),'LineWidth',1)
% t = strcat('Third Level (',num2str(altitudes(3)),'-',num2str(1000),'km)');
% title(t);
% xlabel('Time (days)');
% ylabel('Total Debris');

end
%%
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


