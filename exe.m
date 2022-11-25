%% user inputs
clc; clear; close all;
seasonalAvailability = 1; % 0 means heating and cooling available year-round 
sOaMin = 10; % percent minimum outdoor air damper position setpoint
dataDirectory = 'C:\Users\burak\Dropbox\Research\Maintenance\FDD\hard and soft\buildSys 2022\data\2018';
start = 8;
stop = 20;
wkndOp = 1;
htgSp = 22;
clgSp = 23.5;
htgSpUnocc = 18;
clgSpUnocc = 27;

[energy, znMdlPrmtr] = faultImpactAnalyzer(seasonalAvailability,sOaMin,dataDirectory,start,stop,wkndOp,htgSp,clgSp,htgSpUnocc,clgSpUnocc);
