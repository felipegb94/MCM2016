numlevels = 2;
%given parameters/rate constants
klaunch = 0.0*ones(numlevels,1);   
kcolcraftfrag = 0*2*10^-9;
kcolfragfrag = 0*1*10^-9;
kdecay = .0005;

decaytimes = round(10.^(3*(1:numlevels)./numlevels))';

% number of days to simulate
tend = 200;
% number of timesteps to preallocate
tsteps = 100*tend;

% reaction vectors in columns
% -spacecraft launch
% -catastrophic collision
% -fragment collision
% -orbital decay
%
% the rows are are
% -spacecraft
% -debris
zeta = [1  -1   0   0;
        10 100 100 0];
rxnum = size(zeta,2);
% preallocate for loop
% X tracks population of each altitude level over time
% in the format
% Level 1 spacecraft
% Level 1 debris
% Level 2 spacecraft
% Level 2 debris
% ...

% initial number of each object
X0 = zeros(2*numlevels,1);
% spacecraft
X0(1:2:2*numlevels) = 0;
%fragments
X0(2:2:2*numlevels) = 5000;

% preallocate tracking number of objects
X = zeros(2*numlevels,tsteps);
X(:,1) = X0;
% seperate set of intensity functions for each level
lam = zeros(numlevels,rxnum);
% create birth tracking
pret = max(decaytimes); %extra time before initialization
birth = zeros(numlevels,10*tend+pret);
birth(:,1:pret) = randi(20,numlevels,pret);
life = zeros(numlevels,1);
% random number tracker
rcount = 1;
% generate random number pool
rs = rand(numlevels,3*tsteps);

%generate starting uniform randoms
Pk = log(1./rand(numlevels,rxnum));
Tk = zeros(numlevels,rxnum);
t = zeros(numlevels,tsteps);
count = 1;
newday = true;
while min(t(:,count))<=tend
    % object decay
    % the first time that a new day is achieved,
    % remove objects that appeared that length ago from that
    % level and move them to the next level down
    if newday
        binds = sub2ind(size(birth),(1:numlevels)',...
                        life-decaytimes+1+pret);
        decaynum = birth(binds);
        % subtract from level
        X(2:2:end,count) = X(2:2:end,count) - decaynum;
        % add to next level
        X(2:2:end,count) = X(2:2:end,count) + [decaynum(2:end);0];
    end
    % stochastic aspect
    % calculate intensity functions
    lam(:,1) = klaunch;
    lam(:,2) = kcolcraftfrag*X(2:2:end,count).*X(1:2:end,count);
    lam(:,3) = kcolfragfrag*(X(2:2:end,count)).^2;
    lam(:,4) = kdecay*X(2:2:end,count);
    % find scaled time
    deltk = (Pk-Tk)./lam;
    % find next reaction
    [del,mu] = min(deltk,[],2);
    % add reaction vector
    for z = 1:numlevels
        X(2*z-1:2*z,count+1) = X(2*z-1:2*z,count) + zeta(:,mu(z));
    end
     
    %compute number of new fragments for birth decay process
    diff = X(2:2:end,count+1)-X(2:2:end,count)+decaynum;
    curbinds = sub2ind(size(birth),(1:numlevels)',life+pret+1);
    birth(curbinds) = birth(curbinds)+[diff(2:end);0];
   
    % update scale times
    Tk = Tk + bsxfun(@times,del,lam);
    % update the triggered reaction P
    Pinds = sub2ind([numlevels,rxnum],(1:numlevels)',mu);
    Pk(Pinds) = Pk(Pinds) + log(1./rs(:,rcount));
    rcount = rcount + 1;
    % update overall time
    t(:,count + 1) = t(:,count) + del;
    count = count + 1;
    % update day ticker for each level
    life = floor(t(:,count));
end


%plot the trajectory
figure(1);
% suptitle('Single Trajectory of Collision/Spacecraft Model')
for ind = 1:numlevels
subplot(numlevels,1,ind)
plotyy(t(ind,1:count),X(2*ind-1,1:count),t(ind,1:count),X(2*ind,1:count));
titl = sprintf('Level %d',ind);
title(titl)
% legend('Spacecraft','Debris')
xlabel('time')
ylabel('Number')
end