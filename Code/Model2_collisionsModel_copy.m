close all;
numlevels = 3;
tend = 6000;
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
fragdebris = 2000;
fraginds = sub2ind(size(rxns),2:2:2*numlevels,1:numlevels);
rxns(fraginds) = fragdebris;
%other reactions

% initialize randoms
rs = rand(1,10*tend);
Pk = log(1./rand(size(rxns,2),1));
Tk = zeros(size(rxns,2),1);
currTimeBin = 1;
currTimeBin2 = currTimeBin;
while currTime < tend-1
    disp(currTime)
    % call stochastic reactions function
    rnd = rs(counter);
    [change,delt,Pk,Tk] =...
        stochevent(X(:,currTimeBin),rxns,Pk,Tk,rnd);
    
    history(:,currIndex) = history(:,currIndex)+change(2:2:end);
    %calculate decay
    [history, currIndex, currDecayIndex,decay,flag] = decaylevels(history, currIndex, currDecayIndex, currTime, delt);
    if delt>=1
        delt=.999;
    end
    currTime = currTime + delt;
%     currTimeBin = currDecayIndex(numlevels);
    currTimeBin = ceil(currTime);
    
    if currTimeBin > (currTimeBin2 + 1)
        disp('a');
        disp(currTimeBin);
        disp('b');
        disp(currTimeBin2);
        error('Error for currTimeBin');
    end
    currTimeBin2 = currTimeBin;
    
    if flag
        X(:,currTimeBin) = X(:,currTimeBin-1) + decay;
    end
    X(:,currTimeBin) = X(:,currTimeBin) + change;
    newX = X(:,currTimeBin)
    population(:,counter) = X(:,currTimeBin);
    time(counter) = currTime;
    counter = counter + 1;
end

figure;
subplot(3,1,1);
semilogy(time(1:counter-1), population(2,1:counter-1), 'LineWidth',1)
t = strcat('First Level (',num2str(altitudes(1)),'-',num2str(altitudes(2)),'km)');
title(t);
xlabel('Time (days)');
ylabel('Total Debris');
subplot(3,1,2);
semilogy(time(1:counter-1), population(4,1:counter-1),'LineWidth',1)
t = strcat('Second Level (',num2str(altitudes(2)),'-',num2str(altitudes(3)),'km)');
title(t);
xlabel('Time (days)');
ylabel('Total Debris');
subplot(3,1,3);
plot(time(1:counter-1), population(6,1:counter-1),'LineWidth',1)
t = strcat('Third Level (',num2str(altitudes(3)),'-',num2str(1000),'km)');
title(t);
xlabel('Time (days)');
ylabel('Total Debris');



