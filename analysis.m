clear all
addpath(genpath('.'));
% load dataset
global FileInfo
load('DataSet_1.ppf','-mat');
% do burstsearch
IB = 1;
IT = 10;
%[start,stop] = CUSUM_BurstSearch(MT{1,1},IB,IT,0.05);
%[start,stop] = SlidingTimeWindow_BurstSearch(MT{1,1},5,500,true);
[start,stop] = InterphotonTime_BurstSearch(MT{1,1},15,160);
%[start,stop] = ChangePoint_BurstSearch(MT{1,1},10);
% load true start stops
start_gt = cellfun(@(x) x(1),PhotonNumbers)';
stop_gt = cellfun(@(x) x(end),PhotonNumbers)';

start_dt = zeros(size(start));
stop_dt = zeros(size(start));
for i = 1:numel(start)
    % find the closes start and record difference
    start_dt(i) = min(abs(start(i)-start_gt));
    stop_dt(i) = min(abs(stop(i)-stop_gt));
end

N_photons = stop-start+1;
L = 5;
%%
figure;hold on;
histogram(start_dt(N_photons>=L),-0.5:1:20.5)
histogram(stop_dt(N_photons>=L),-0.5:1:20.5)
% figure;histogram(stop_gt-start_gt+1);
% figure;histogram(N_photons(N_photons>L))
fprintf('%i bursts detected.\n%i bursts above threshold of  %i photons (%.2f %% of total).\n',...
    numel(N_photons),sum(N_photons >= L),L,100*sum(N_photons >= L)/numel(N_photons));