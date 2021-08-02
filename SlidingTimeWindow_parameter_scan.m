clear all
addpath(genpath('.'));
% load dataset
load('DataSet_1.ppf','-mat');
global FileInfo
FileInfo = Info;
Photons = MT{1,1};
% determine ground truth
load('PhotonBurstData.mat');
start_gt = cellfun(@(x) x(1),PhotonNumbers)';
stop_gt = cellfun(@(x) x(end),PhotonNumbers)';

% do burstsearch
% scan over threshold and false positive percentage
[M,T] = meshgrid(1:1:100,100:100:10000);
%[IT,alpha] = meshgrid(10,[1E-6,0.001,0.01:0.01:0.99]);
%[IT,alpha] = meshgrid(1:1:50,1E-3);
dt = zeros(size(M));
dt_UB = zeros(size(M));
dt_LB = zeros(size(M));
N_bursts = zeros(size(M));
for i = 1:size(M,1)
    for j = 1:size(M,2)
        [start,stop] = SlidingTimeWindow_BurstSearch(Photons,M(i,j),T(i,j),true);
        %[start,stop] = SlidingTimeWindow_BurstSearch(MT{1,1},5,1000,true);
        %[start,stop] = InterphotonTime_BurstSearch(MT{1,1},1,260);
        %[start,stop] = ChangePoint_BurstSearch(MT{1,1},10);
        % load true start stops
        
        start_dt = zeros(size(start));
        stop_dt = zeros(size(start));
        for k = 1:numel(start)
            % find the closes start and record difference
            start_dt(k) = min(abs(start(k)-start_gt));
            stop_dt(k) = min(abs(stop(k)-stop_gt));
        end
        N_bursts(i,j) = numel(start);
        N_photons = stop-start+1;
        L = 10;

        dt(i,j) = mean([start_dt(N_photons>=L); stop_dt(N_photons>=L)]);
        dt_UB(i,j) = prctile([start_dt(N_photons>=L); stop_dt(N_photons>=L)],1-0.32/2);
        dt_LB(i,j) = prctile([start_dt(N_photons>=L); stop_dt(N_photons>=L)],0.32/2);
    end
    disp(i*size(M,2)/numel(M));
end
%%
figure;
contourf(M,T,dt,'LevelList',0:1:10);
%surf(M./T*1000,T,N_bursts);
%loglog(M,dt,'LineWidth',1);
%set(gca,'YScale','log');