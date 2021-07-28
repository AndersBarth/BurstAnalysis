clear all
addpath(genpath('.'));
% load dataset
load('DataSet_1.ppf','-mat');
global FileInfo
FileInfo = Info;
% do burstsearch
IB = 1;
IT = 10;
%[start,stop] = CUSUM_burstsearch(MT{1,1},IB,IT);
%[start,stop] = SlidingTimeWindow_BurstSearch(MT{1,1},5,500,true);
%[start,stop] = InterphotonTime_BurstSearch(MT{1,1},10,100);
[start,stop] = ChangePoint_BurstSearch(MT{1,1},10);
% load true start stops
load('PhotonBurstData.mat');
start_gt = cellfun(@(x) x(1),PhotonNumbers)';
stop_gt = cellfun(@(x) x(end),PhotonNumbers)';

start_dt = zeros(size(start));
stop_dt = zeros(size(start));
for i = 1:numel(start)
    % find the closes start and record difference
    start_dt(i) = min(abs(start(i)-start_gt));
    stop_dt(i) = min(abs(stop(i)-stop_gt));
end
figure;histogram(start_dt,0:1:20);hold on;histogram(stop_dt,0:1:20)
disp(mean(start_dt)-mean(stop_dt));