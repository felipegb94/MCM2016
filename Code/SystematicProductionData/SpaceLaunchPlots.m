clear;
clf;
close all;

%% Get Launch Target height
% Year, number of flights.
load('LEOLaunchData.mat');
% Year, unmaned launches, manned, total
load('totalLaunchData.mat');

time = LEOLaunchData(:,1);
LEOLaunches = LEOLaunchData(:,2);
totalLaunches = totalLaunchData(:,4);

figure;
hold on;
grid on;
title('Space Launches Per Year')
p = plot(time, totalLaunches, time, LEOLaunches);
p(1).Marker = '*';
p(2).Marker = '*';

xlabel('Year');
ylabel('Launches');


load('LEOData.mat');
ISS = LEOData.ISS;
S = LEOData.S;
Other = LEOData.Other;

SumISS = sum(ISS);
SumS = sum(S);
SumOther = sum(Other);
SumLEOLaunches = sum(LEOLaunches(11:16));

RatioISSLaunches = SumISS(1,2) / SumLEOLaunches
RatioSLaunches = SumS(1,2) / SumLEOLaunches
RatioOtherLaunches = SumOther(1,2) / SumLEOLaunches

AverageLEOLaunches = mean(LEOLaunches)

% We want the ratios to add up to 1. Currently they add ~98.5
sumRatios = RatioISSLaunches + RatioSLaunches + RatioOtherLaunches;
toAdd = (1 - sumRatios)/3;

ratios = [RatioISSLaunches;RatioSLaunches;RatioOtherLaunches]
ratios = ratios + toAdd;
Launches = floor(ratios * AverageLEOLaunches);


% Plot averages
plot(time, repmat(AverageLEOLaunches,size(time,1),1));
plot(time, repmat(AverageLEOLaunches * ratios(1),size(time,1),1));
plot(time, repmat(AverageLEOLaunches * ratios(2),size(time,1),1));
plot(time, repmat(AverageLEOLaunches * ratios(3),size(time,1),1));

legend('Total Launches', 'Total Launches into LEO', 'Average Launches into LEO', 'Average ISS Launches', 'Average Sun Synchronous LEO', 'Average Other LEO Launches');
axis([2000,2015,0,100]);

LaunchData.AverageLEOLaunches = AverageLEOLaunches;
LaunchData.ratios = ratios;
LaunchData.ISSLaunches = floor(AverageLEOLaunches * ratios(1));
LaunchData.SLaunches = floor(AverageLEOLaunches * ratios(2));
LaunchData.OtherLaunches = floor(AverageLEOLaunches * ratios(3));

save('../data/LaunchData.mat','LaunchData');

%% Get Launch heights
rng('shuffle')
% 
HeightRangeISS = [370, 460]; 
LaunchHeightsISS = HeightRangeISS(1) + rand(1, Launches(1)) * (HeightRangeISS(2) - HeightRangeISS(1));

HeightRangeS = [500, 800];
LaunchHeightsS = HeightRangeS(1) + rand(1, Launches(2)) * (HeightRangeS(2) - HeightRangeS(1));

HeightRangeOthers = [200, 1000];
LaunchHeightsOthers = HeightRangeOthers(1) + rand(1, Launches(3)) * (HeightRangeOthers(2) - HeightRangeOthers(1));

LaunchHeights = [LaunchHeightsISS, LaunchHeightsS, LaunchHeightsOthers]
figure;
hist(LaunchHeights)
xlabel('Altitude (km)');
ylabel('Num Satellites (km)');

figure;
load('../data/satelliteAltitudes.mat')
hist(satelliteAltitudes)
xlabel('Altitude (km)');
ylabel('Num Satellites (km)');











