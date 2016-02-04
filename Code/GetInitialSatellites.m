function [initialSatellites ] = GetInitialSatellites(numlevels)
%GetInitialOperationalSatellites 

load('data/satelliteAltitudes.mat');
initialSatellites = hist(satelliteAltitudes,numlevels);

end

