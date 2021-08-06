function SlidingTimeWindow_parameter_scan(filename,L)
% L - Minimum number of photons per burst
if nargin < 2
    L = 5;
end

% scan over threshold and false positive percentage
M_range = 2:1:25;
T_range = 100:100:2000;
[M,T] = meshgrid(M_range,T_range);

% load dataset
load(filename,'-mat');

% determine ground truth
start_gt = cellfun(@(x) x(1),BurstPhotonNumbers)';
stop_gt = cellfun(@(x) x(end),BurstPhotonNumbers)';
N_photons_gt = stop_gt - start_gt + 1;

% do burstsearch
dt = zeros(size(M));
dt_UB = zeros(size(M));
dt_LB = zeros(size(M));
N_bursts = zeros(size(M));
N_FP = zeros(size(M));
N_P = zeros(size(M));
N_TP = zeros(size(M));
N_split = zeros(size(M));
for i = 1:size(M,1)
    for j = 1:size(M,2)
        [start,stop] = SlidingTimeWindow_BurstSearch(MT,M(i,j),T(i,j),true);
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
        % number of detected bursts
        N_bursts(i,j) = numel(start);
        N_photons = stop-start+1;
        
        % mean photon deviation
        L = 5;
        dt(i,j) = mean([start_dt(N_photons>=L); stop_dt(N_photons>=L)]);
        %dt_UB(i,j) = prctile([start_dt(N_photons>=L); stop_dt(N_photons>=L)],1-0.32/2);
        %dt_LB(i,j) = prctile([start_dt(N_photons>=L); stop_dt(N_photons>=L)],0.32/2);
        
        % filter based on number of photons
        start = start(N_photons>=L);
        stop = stop(N_photons>=L);
        % identify false positives outside of the burst regions
        FP = false(size(start));
        for k = 1:numel(start)
            % find index of last ground-truth stop before start
            ix_stop = find(stop_gt < start(k),1,'last');
            if isempty(ix_stop)
                ix_stop = 0;
            end
            % find index of next ground-truth start after stop
            ix_start = find(start_gt > stop(k),1,'first');
            if isempty(ix_start)
                ix_start = numel(start_gt)+1;
            end
            if ix_start-ix_stop == 1
                % this means the burst starts after the stop of
                % ground-truth burst i and ends before the start of
                % ground-truth burst i+1
                FP(k) = true;
            end
        end
        N_FP(i,j) = sum(FP); % number of false positives
        N_P(i,j) = numel(start); % total number of detected bursts
        
        % find number of bursts per true burst region (true positives and splitting)
        n_per_burst = zeros(size(start_gt));
        for k = 1:numel(start_gt)
            n_per_burst(k) = sum(stop >= start_gt(k) & start <= stop_gt(k));
        end        
        % find true positives, i.e. true bursts that were detected by at least one burst
        N_TP(i,j) = sum(n_per_burst > 0);
        N_split(i,j) = sum(n_per_burst)./N_TP(i,j); % average number of bursts per detected true positive
    end
    disp(i*size(M,2)/numel(M));
end

% Quantities related to the sensitivity and specificity
%       We do know:
% P   - The number of real bursts in the data. P = TP + FN
% TP  - True positives, i.e. the number of "hits" where we found at least one burst. TP = P - FN
% FP  - False positives, i.e. bursts detected outside of true regions, "false alarm". FP = N - TN
% FN  - False negatives, i.e. true bursts that where not detected. FN = P - TP
% These quantities are not defined here as we do not simulate "negative events":
% N   - The number of negatives in the data. N = TN + FP
% TN  - True negatives, i.e. the not detecting a burst outside of burst regions.
%
%       From the known quantities, we can estimate:
% TPR - True positive rate (sensitivity, hit rate). TPR = TP/P = 1-FNR
% FNR - False negative rate (miss rate). FNR = FN/P = 1-TPR
% FDR - False discovery rate. FDR = FP/FP+TP = 1-PPV
%       Inverse quantities to these are:
% PPV - Positive predictive value (precision). PPV = TP/TP+FP = 1-FDR
%       We can not estimate these quantities:
% TNR - True negative rate (specificity, selectivity). TNR = TN/N
% NPV - Negative predictive value. NPV = TN/TN+FN
% FPR - False positive rate (fall-out). FPR = FP/N
% FOR - False omission rate. FOR = FN/FN+TN

% the number of real bursts in the data
P = numel(start_gt);
% true positive rate (sensitivity, hit rate)
TPR = N_TP./P;
% false negative rate (miss rate)
FNR = 1 - TPR;
% false discovery rate
FDR = N_FP./(N_FP+N_TP);

%% plot and save figure
lw = 1;
fs = 8;
figure('Color',[1,1,1],'Units','inch','Position',[0,0,9.5,5]); hold on;
tiledlayout(2,3);

nexttile; hold on; colormap('vik');
imagesc(M(1,:),T(:,1),dt,'AlphaData',dt>0);
contour(M,T,N_bursts,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,5],'YLim',[0,T(end,1)],'Box','on','LineWidth',lw,'FontSize',fs);
xlabel('photons per time window, M');
ylabel('time window (\mus)');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
[~,ix] = min(dt,[],1,'linear');
plot(M(ix),T(ix),'.w');
%plot(M_range,M_range*mean(T(ix)./M(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(M(ix)./T(ix))*1E3));
c = colorbar;
c.Label.String = 'photon deviation';
title('photon deviation');
ax = gca;
text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);

nexttile; hold on;
imagesc(M(1,:),T(:,1),TPR,'AlphaData',~isnan(TPR));
%contour(M,T,TPR,'LevelList',1,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,1],'Box','on','LineWidth',lw,'FontSize',fs);
c = colorbar;
c.Label.String = 'true positive rate';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('true positive rate');
ax = gca;
text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);

nexttile; hold on;
imagesc(M(1,:),T(:,1),FDR,'AlphaData',~isnan(FDR));
%contour(M,T,FDR,'LevelList',0,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,1],'Box','on','LineWidth',lw,'FontSize',fs);
c = colorbar;
c.Label.String = 'false discovery rate';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('false discovery rate');
ax = gca;
text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);

nexttile; hold on;
imagesc(M(1,:),T(:,1),N_P);
contour(M,T,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,2500],'Box','on','LineWidth',lw,'FontSize',fs);
c = colorbar;
c.Label.String = 'number of bursts';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('number of detected bursts');
ax = gca;
text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);

nexttile; hold on;
imagesc(M(1,:),T(:,1),N_split,'AlphaData',~isnan(FDR));
set(gca,'Box','on','LineWidth',lw,'FontSize',fs);
c = colorbar;
c.Label.String = 'splitting frequency';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('splitting frequency');
ax = gca;
text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);


print(gcf,[filename(1:end-3) 'png'],'-dpng','-painters','-r300');
print(gcf,[filename(1:end-3) 'eps'],'-depsc','-painters');