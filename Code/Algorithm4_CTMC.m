%given parameters/rate constants
k1 = 200;   
k2 = 10;
k3 = .01;
dm = 25;
dp = 1;
dd = 1;
% number of iterations- Monte Carlo
n = 1;

% reaction vectors in columns
zeta = [0  0  0  0  0  0;
        1  0 -1  0  0  0;
        0  1  0 -1 -2  0;
        0  0  0  0  1 -1];

% preallocate for loop
X0 = [1 10 50 10]';
X = zeros(4,5000);
XT5 = zeros(4,n);
X(:,1) = X0;
lam = zeros(6,1);
% random number tracker
rcount = 1;
% generate random number pool
rs = rand(1,50000*n);
for idx = 1:n
%generate starting uniform randoms
Pk = log(1./rand(6,1));
Tk = zeros(6,1);
t = zeros(1,5000);
count = 1;
while t(count)<5
    % calculate intensity functions
    lam(1) = k1;
    lam(2) = k2*X(2,count);
    lam(3) = dm*X(2,count);
    lam(4) = dp*X(3,count);
    lam(5) = k3*X(3,count)*(X(3,count)-1);
    lam(6) = dd*X(4,count);
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
% final value to tracker
XT5(:,idx) = X(:,count-1);
end
%truncate
X = X(:,1:count);
t = t(1:count);
%plot the last trajectory
figure(1)
plot(t,X)
title('Single Trajectory of Gene Transcription/Translation')
legend('Gene','mRNA','Protein','Dimer')
xlabel('time')
ylabel('Number of Species')

% calculate expectation of each species
Xmean = mean(XT5,2);

% calculate 2nd moment of each species
Xsecond = mean(XT5.^2,2);
Xvar = Xsecond - Xmean.^2;
display(Xmean(4),'Mean')
display(Xvar(4),'Variance')