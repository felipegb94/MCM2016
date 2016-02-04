function [change,delt,Pk,Tk] = stochevent(X,rxns,Pk,Tk,rnd)
% calculate the next step in the stochastic process

%intensity parameters
kfragfrag = 1*10^-6;
kdummy = 3;


%calculate intensity functions
lam = zeros(size(rxns,2),1);
% fragment-fragment collisions
lam(1:end-1) = kfragfrag.*X(2:2:end);
%dummy intensity function
lam(end) = kdummy;

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