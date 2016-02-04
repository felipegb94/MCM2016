numlevels = 50;
%given parameters/rate constants
klaunch = 0.1*ones(numlevels,1);   
kcolcraftfrag = 2*10^-8;
kcolfragfrag = 1*10^-9;
kdecay = 3*ones(numlevels,1);

% number of timesteps
tsteps = 200;

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
        10 100 100 -2];
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
X0 = repmat([200 5000]',numlevels,1);
% preallocate tracking number of objects
X = zeros(2*numlevels,tsteps);
X(:,1) = X0;
% seperate set of intensity functions for each level
lam = zeros(numlevels,rxnum);

% random number tracker
rcount = 1;
% generate random number pool
rs = rand(numlevels,10000);

%generate starting uniform randoms
Pk = log(1./rand(numlevels,rxnum));
Tk = zeros(numlevels,rxnum);
t = zeros(numlevels,tsteps);
count = 1;
while count<=tsteps
    % calculate intensity functions
    lam(:,1) = klaunch;
    lam(:,2) = kcolcraftfrag*X(2:2:end,count).*X(1:2:end,count);
    lam(:,3) = kcolfragfrag*(X(2:2:end,count)).^2;
    lam(:,4) = kdecay;
    % find scaled time
    deltk = (Pk-Tk)./lam;
    % find next reaction
    [del,mu] = min(deltk,[],2);
    % add reaction vector
    for z = 1:numlevels
        X(2*z-1:2*z,count+1) = X(2*z-1:2*z,count) + zeta(:,mu(z));
    end
    % update scale times
    Tk = Tk + bsxfun(@times,del,lam);
    % update the triggered reaction P
    Pinds = sub2ind([numlevels,rxnum],(1:numlevels)',mu);
    Pk(Pinds) = Pk(Pinds) + log(1./rs(:,rcount));
    rcount = rcount + 1;
    % update overall time
    t(:,count + 1) = t(:,count) + del;
    count = count + 1; 
end


%plot the trajectory
figure;
for ind = 1:numlevels
subplot(numlevels,1,ind)
plotyy(t(ind,:),X(2*ind-1,:),t(ind,:),X(2*ind,:));
titl = sprintf('Level %d',ind);
title(titl)
% legend('Spacecraft','Debris')
end



xlabel('time')
ylabel('Number')
suptitle('Single Trajectory of Collision/Spacecraft Model')

plotyy(t(ind,:),X(2*ind-1,:),t(ind,:),X(2*ind,:));
titl = sprintf('Level %d',ind);
title(titl)
%plotyy(t(4,:),X(7,:),t(4,:),X(8,:));

