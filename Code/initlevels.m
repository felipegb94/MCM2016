function [history, decaytimes, offset,altitudes, collisionProbs] = initlevels(numlevels,tdays)

basealt = 200;
maxalt = 1600;
fine = 1000;
earthR = 6371;

%calculate decay times
load('data/altdecaytimes_extrapolated')
altitudes = round(linspace(basealt,maxalt,numlevels+1))';
decaytimes = ceil(interp1(altdecaytimes_extrapolated(:,1),...
                altdecaytimes_extrapolated(:,2),altitudes));
decaytimes = decaytimes(1:end-1);

%determine density and quantity values
load('data/altdensity')
altdensity(:,2) = 13*altdensity(:,2);
altitudes2 = round(linspace(basealt,maxalt,fine*numlevels))';
densities = interp1(altdensity(:,1),altdensity(:,2),...
                    altitudes2);
                
numdebris = zeros(numlevels,1);


for d_ind = 1:length(densities)-1
   vol=4/3*pi*((altitudes2(d_ind+1)+earthR)^3 - (altitudes2(d_ind)+earthR)^3);
   temp = floor((d_ind+fine)/fine);
   numdebris(temp) =numdebris(temp)+vol*densities(d_ind+1);
end
numdebris = ceil(numdebris);




%temp values
% decaytimes = round(linspace(1,1000,numlevels))';

% create initial array
history = zeros(numlevels,max(decaytimes)+tdays+2);
offset = max(decaytimes);
for ind = 1:numlevels
    debrisdistr = randi([offset-decaytimes(ind)+1,offset],1,numdebris(ind));
    for ind2 = 1:length(debrisdistr)
        history(ind,debrisdistr(ind2)) = history(ind,debrisdistr(ind2))+1;
    end
end




end