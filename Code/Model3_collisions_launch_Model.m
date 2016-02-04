function Model3_collisions_launch_Model

numlevels = 6;
tend = 36500;
currTime = 0;
[history, decaytimes, offset, altitudes] = initlevels(numlevels,tend);

currDecayIndex = offset + 1 - decaytimes;
currIndex = offset + 1;

counter = 1;
population = zeros(2*numlevels,5*tend);
time = zeros(1,5*tend);


% prepare stochastic variables
X = zeros(2*numlevels,tend);
X(1:2:end,1) = GetInitialSatellites(numlevels);
X(2:2:end,1) = sum(history,2);
% create reaction vectors
rxns = zeros(2*numlevels,4*numlevels+1);
% fragment-fragment reactions
fragdebris = 20;
fraginds = sub2ind(size(rxns),2:2:2*numlevels,1:numlevels);
rxns(fraginds) = fragdebris;
% fragment-large satellite reactons
spfragdebris = 10000;
spfraginds = sub2ind(size(rxns),2:2:2*numlevels,numlevels+1:2*numlevels);
rxns(spfraginds) = spfragdebris;
spcrashinds = sub2ind(size(rxns),1:2:2*numlevels,numlevels+1:2*numlevels);
rxns(spcrashinds) = -1;
% launches
launchinds = sub2ind(size(rxns),1:2:2*numlevels,2*numlevels+1:3*numlevels);
rxns(launchinds) = 1;
launchdebinds = sub2ind(size(rxns),2:2:2*numlevels,2*numlevels+1:3*numlevels);
rxns(launchdebinds) = 70;
LaunchNums = GetYearlyLaunchDebris(numlevels);
% fragment removal
% calculate number of fragments to remove
remlevels = [];
for indx = 1:numlevels
   if altitudes(indx+1)>600 && altitudes(indx)<1200
        remlevels = [remlevels indx];
   end
end
removalinds = sub2ind(size(rxns),2*remlevels,3*numlevels+remlevels);
rxns(removalinds) = -round(356/length(remlevels));

%calculate volumes for density scaling
vol = zeros(numlevels,1);
for lv = 1:numlevels
   vol(lv) = 4/3*pi*((altitudes(lv+1)+6371)^3 - (altitudes(lv)+6371)^3);
end
% initialize randoms
rs = rand(1,1000*tend);
Pk = log(1./rand(size(rxns,2),1));
Tk = zeros(size(rxns,2),1);
currTimeBin = 1;
collcount = 0;
fragcollcount = 0;
collt = [];
while currTime < tend-1
%     if(mod(floor(currTime),1000) == 0)
%         disp(currTime)
%     end
    if mod(floor(currTime),365)==0
        LaunchNums = GetYearlyLaunchDebris(numlevels);
    end
    %cutoff removal
    if currTimeBin == 1825
       rxns(rxns==-89) = 0; 
    end
    % call stochastic reactions function
    rnd = rs(counter);
    [change,delt,Pk,Tk,rxnum] =...
        stochevent_model3(X(:,currTimeBin),rxns,Pk,Tk,rnd,vol,LaunchNums);
    if rxnum > numlevels && rxnum <= 2*numlevels
        collcount = collcount+1;
        collind = rxnum-numlevels;
        colltrack = [currTimeBin;collind];
        collt = [collt colltrack];
    elseif rxnum <= numlevels
        fragcollcount = fragcollcount + 1;
    end
    
    if rxnum <= 3*numlevels
        history(:,currIndex) = history(:,currIndex)+change(2:2:end);
    end
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
    % correct for negatives
    
    X(X(:,currTimeBin)<0,currTimeBin) = 0;
%     newX = X(:,currTimeBin);
    population(:,counter) = X(:,currTimeBin);
    time(counter) = currTime;
    counter = counter + 1;
end

% Plotting %
%%%%%%%%%%%%
figure(2);
for plotind = 1:numlevels
subplot(numlevels,2,2*plotind-1);
plot(time(1:counter-1), population(2*plotind,1:counter-1), 'LineWidth',1)
t = sprintf('Level %d: %d - %d Km',plotind,altitudes(plotind),altitudes(plotind+1));
title(t);
xlabel('Time (days)');
ylabel('Total Debris');
if ~isempty(collt)
    inds = find(collt(2,:)==plotind);
    hold on
    plot(collt(1,inds),zeros(1,length(inds)),'r*')
    hold off
end

subplot(numlevels,2,2*plotind);
plot(time(1:counter-1), population(2*plotind-1,1:counter-1), 'LineWidth',1)
t = sprintf('Level %d: %d - %d Km',plotind,altitudes(plotind),altitudes(plotind+1));
title(t);
xlabel('Time (days)');
ylabel('Total Spacecraft');
if ~isempty(collt)
    inds = find(collt(2,:)==plotind);
    hold on
    plot(collt(1,inds),zeros(1,length(inds)),'r*')
    hold off
end
end
% function return
removalamount = -2*sum(sum(rxns(1:2*numlevels,3*numlevels+1:4*numlevels)))
initials = X(:,1);
finals = X(:,currTimeBin);
totinitials = [sum(initials(1:2:end)) sum(initials(2:2:end))]
totfinals = [sum(finals(1:2:end)) sum(finals(2:2:end))]
collcount = [collcount fragcollcount]
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