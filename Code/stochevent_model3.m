function [change,delt,Pk,Tk,mu] = ...
    stochevent_model3(X,rxns,Pk,Tk,rnd,vols,LaunchNums)
% calculate the next step in the stochastic process
numlevels = length(X)/2;
%calculate density scaling factors
%intensity parameters
collprobs = GetNumCollisions(X,vols);
klaunches = LaunchNums./365;
kremoval = 2/numlevels;
kdummy = 1;



%calculate intensity functions
lam = zeros(size(rxns,2),1);
% fragment-fragment collisions
lam(1:numlevels) = collprobs(:,1);
% spacecraft-fragment collisions
lam(numlevels+1:2*numlevels) = collprobs(:,2);
% spacecraft launches
lam(2*numlevels+1:3*numlevels) = klaunches;
% fragment removal
lam(3*numlevels+1:4*numlevels) = kremoval;
%dummy intensity function
lam(end) = kdummy;
% correct for negative lambdas
lam(lam<0) = 0;

% find scaled time
deltk = (Pk-Tk)./lam;
% find next reaction
[delt,mu] = min(deltk,[],1);
%select reaction vector and return to call
change = rxns(:,mu);

%update scale times
Tk = Tk + delt*lam;
% update the triggered reaction P
Pk(mu) = Pk(mu) + log(1/rnd);
end