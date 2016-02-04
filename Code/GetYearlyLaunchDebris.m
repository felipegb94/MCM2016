function LaunchNumbers = GetYearlyLaunchDebris( numlevels )
%GetYearlyLaunchDebris 
%   Compute the amount of mission-related debris on a year. The
%   computations in this function are data-driven, this means that the
%   average number of launches and height estimates were obtained from
%   actual data.

rng('shuffle')
load('data/LaunchData.mat')
LaunchesISS = LaunchData.ISSLaunches;
LaunchesS = LaunchData.SLaunches;
LaunchesOther = LaunchData.OtherLaunches;

% 
HeightRangeISS = [370, 460]; % ISS height range 
LaunchHeightsISS = HeightRangeISS(1) + rand(1, LaunchesISS) * (HeightRangeISS(2) - HeightRangeISS(1));

HeightRangeS = [600, 800]; % Sun synchronous orbit average heights
LaunchHeightsS = HeightRangeS(1) + rand(1, LaunchesS) * (HeightRangeS(2) - HeightRangeS(1));

HeightRangeOthers = [200, 1600]; % Other uncategorized launches.
LaunchHeightsOthers = HeightRangeOthers(1) + rand(1, LaunchesOther) * (HeightRangeOthers(2) - HeightRangeOthers(1));

LaunchHeights = [LaunchHeightsISS, LaunchHeightsS, LaunchHeightsOthers];

LaunchNumbers = hist(LaunchHeights,numlevels)



end

