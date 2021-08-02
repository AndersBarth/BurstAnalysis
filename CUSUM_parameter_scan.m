%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For the CUSUM, the interesting parameter seems to be the error rate %%%
%%% alpha. By default, it is set to 1/Nphotons, but the resulting low   %%%
%%% rate causes a high detection delay.                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
IB = 1;
% scan over threshold and false positive percentage
[IT,alpha] = meshgrid(10,logspace(-9,-0.05,100));
%[IT,alpha] = meshgrid(10,[1E-6,0.001,0.01:0.01:0.99]);
%[IT,alpha] = meshgrid(1:1:50,1E-3);
dt = zeros(size(IT));
dt_UB = zeros(size(IT));
dt_LB = zeros(size(IT));
N_bursts = zeros(size(IT));
for i = 1:size(IT,1)
    for j = 1:size(IT,2)
        [start,stop] = CUSUM_BurstSearch(Photons,IB,IT(i,j),alpha(i,j));
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
    disp(i*size(IT,2)/numel(IT));
end
%%
figure;
loglog(alpha,dt,'LineWidth',1);
%set(gca,'YScale','log');