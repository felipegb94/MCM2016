%given parameters/rate constants
klaunch = 0;   
kcolcraftfrag = 2*10^-8;
kcolfragfrag = 1*10^-9;
kdecay = 5;

% number of timesteps
tsteps = 1000;

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

% preallocate for loop
X0 = [200 5000]';
X = zeros(2,tsteps);
X(:,1) = X0;
lam = zeros(4,1);
% random number tracker
rcount = 1;
% generate random number pool
rs = rand(1,10000);

%generate starting uniform randoms
Pk = log(1./rand(4,1));
Tk = zeros(4,1);
t = zeros(1,tsteps);
count = 1;
while t(end)==0
    % calculate intensity functions
    lam(1) = klaunch;
    lam(2) = kcolcraftfrag*X(2,count)*X(1,count);
    lam(3) = kcolfragfrag*(X(2,count))^2;
    lam(4) = kdecay;
    % find scaled time
    deltk = (Pk-Tk)./lam;
    % find next reaction
    [del,mu] = min(deltk);
    % add reaction vector
    X(:,count+1) = X(:,count) + zeta(:,mu);
    % update scale times
    Tk = Tk + del*lam;
    % update the triggered reaction P
    Pk(mu) = Pk(mu) + log(1/rs(rcount));
    rcount = rcount + 1;
    % update overall time
    t(count + 1) = t(count) + del;
    count = count + 1; 
end


%plot the trajectory
figure(1)
plotyy(t,X(1,:),t,X(2,:))
title('Single Trajectory of Collision/Spacecraft Model')
legend('Spacecraft','Debris')
xlabel('time')
ylabel('Number')