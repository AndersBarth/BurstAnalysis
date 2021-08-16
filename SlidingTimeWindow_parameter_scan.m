function SlidingTimeWindow_parameter_scan(foldername,L)
% L - Minimum number of photons per burst
if nargin < 2
    L = 5;
end
time_res = 12.5E-9;% scan over threshold and false positive percentage
M_range = 2:1:25;
T_range = 100:100:2000;
[M,T] = meshgrid(M_range,T_range);
T_time_res = T*1e-6/time_res;

% initialize averaged variables over all files
N_P_total = zeros(size(M)); % number of detected bursts with at least L photons
dt_total = zeros(size(M)); % mean photon deviation of start/stop
dt_start_total = zeros(size(M)); % mean photon deviation of start
dt_stop_total = zeros(size(M)); % mean photon deviation of stop
dt_start_array = cell(size(M)); % array of all photon deviations of starts
dt_stop_array = cell(size(M)); % array of all photon deviations of stops
dt_time_total = zeros(size(M)); % mean time deviation of start/stop
dt_time_start_total = zeros(size(M)); % mean time deviation of start
dt_time_stop_total = zeros(size(M)); % mean time deviation of stop
dt_time_start_array = cell(size(M)); % array of all time deviations of starts
dt_time_stop_array = cell(size(M)); % array of all time deviations of stops
N_FP_total = zeros(size(M)); % false positives
N_TP_total = zeros(size(M)); % true positives
N_split_total = zeros(size(M)); % splitting frequency
P_total = 0; % total number of ground-truth positives (i.e. bursts)

files = dir([foldername filesep '*.ppf']);
%files = files(1); % for testing, only take one file
fprintf('Analyzing %i files in current folder...\n',numel(files));
ll = fprintf('0 %%');
for f = 1:numel(files)
    % load dataset
    load(fullfile(files(f).folder,files(f).name),'-mat');
    % filter empty bursts (true bursts that did not contain any photons)
    BurstPhotonNumbers = BurstPhotonNumbers(cellfun(@(x) ~isempty(x),BurstPhotonNumbers));
    % determine ground truth
    start_gt = cellfun(@(x) x(1),BurstPhotonNumbers)';
    stop_gt = cellfun(@(x) x(end),BurstPhotonNumbers)';
    N_photons_gt = stop_gt - start_gt + 1;

    % do burstsearch
    dt = zeros(size(M));
    dt_start = zeros(size(M));
    dt_stop = zeros(size(M));
    dt_time = zeros(size(M));
    dt_start_time = zeros(size(M));
    dt_stop_time = zeros(size(M));
    
    N_FP = zeros(size(M));
    N_P = zeros(size(M));
    N_TP = zeros(size(M));
    N_split = zeros(size(M));
    for i = 1:size(M,1)
        for j = 1:size(M,2)
            [start,stop] = SlidingTimeWindow_BurstSearch(MT,M(i,j),T_time_res(i,j),true);
            %[start,stop] = SlidingTimeWindow_BurstSearch(MT{1,1},5,1000,true);
            %[start,stop] = InterphotonTime_BurstSearch(MT{1,1},1,260);
            %[start,stop] = ChangePoint_BurstSearch(MT{1,1},10);
            
            % number of photons per burst
            N_photons = stop-start+1;
            % filter based on number of photons
            start = start(N_photons>=L);
            stop = stop(N_photons>=L);
            
             % load true start stops
            start_dt = zeros(size(start));
            stop_dt = zeros(size(start));
            start_dt_time = zeros(size(start));
            stop_dt_time = zeros(size(start));
            for k = 1:numel(start)
                % find the closest start/stop
                [~, start_dt_ix] = min(abs(start(k)-start_gt));
                [~, stop_dt_ix] = min(abs(stop(k)-stop_gt));
                % record photon difference
                start_dt(k) = start(k)-start_gt(start_dt_ix);
                stop_dt(k) = stop(k)-stop_gt(stop_dt_ix);
                % record time difference
                start_dt_time(k) = MT(start(k))-MT(start_gt(start_dt_ix));
                stop_dt_time(k) = MT(stop(k))-MT(stop_gt(stop_dt_ix));
            end
            % mean photon deviation
            dt_start(i,j) = mean(start_dt);
            dt_stop(i,j) = mean(stop_dt);
            dt_start(isnan(dt_start)) = 0;
            dt_stop(isnan(dt_stop)) = 0;
            dt(i,j) = sqrt(mean([start_dt;stop_dt].^2));
            dt_start_array{i,j} = [dt_start_array{i,j}; start_dt];
            dt_stop_array{i,j} = [dt_stop_array{i,j}; stop_dt];
            % mean time deviation
            dt_start_time(i,j) = mean(start_dt_time);
            dt_stop_time(i,j) = mean(stop_dt_time);
            dt_start_time(isnan(dt_start_time)) = 0;
            dt_stop_time(isnan(dt_stop_time)) = 0;
            dt_time(i,j)= sqrt(mean([start_dt_time;stop_dt_time].^2));
            dt_time_start_array{i,j} = [dt_time_start_array{i,j}; start_dt_time];
            dt_time_stop_array{i,j} = [dt_time_stop_array{i,j}; stop_dt_time];
            
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
            N_split(i,j) = sum(n_per_burst); % number of bursts within true positive regions
        end
        %disp(i*size(M,2)/numel(M));
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

    % sum up quantities
    N_P_total = N_P_total + N_P; % number of detected bursts with at least L photons
    dt_total = dt_total + dt; % mean photon deviation of start/stop
    dt_start_total = dt_start_total + dt_start; % mean photon deviation of start
    dt_stop_total = dt_stop_total + dt_stop; % mean photon deviation of stop
    dt_time_total = dt_time_total + dt_time; % mean time deviation of start/stop
    dt_time_start_total = dt_time_start_total + dt_start_time; % mean time deviation of start
    dt_time_stop_total = dt_time_stop_total + dt_stop_time; % mean time deviation of stop
    N_FP_total = N_FP_total + N_FP; % number of false positives
    N_TP_total = N_TP_total + N_TP; % number of true positives
    N_split_total = N_split_total + N_split; % splitting frequency
    P_total = P_total + P; % total number of true bursts
    
    fprintf(repmat('\b',1,ll));
    ll = fprintf('%i %%',round(100*f/numel(files)));
end
n_files = numel(files);
dt = dt_total/n_files;
dt_start = dt_start_total/n_files;
dt_stop = dt_stop_total/n_files;
dt_time = dt_time_total/n_files;
dt_time_start = dt_time_start_total/n_files;
dt_time_stop = dt_time_stop_total/n_files;
N_split = N_split_total./N_TP_total;
N_P = N_P_total/n_files;
P = P_total/n_files;
 % true positive rate (sensitivity, hit rate)
TPR = N_TP_total./P_total;
% false negative rate (miss rate)
% FNR = 1 - TPR;
% false discovery rate
FDR = N_FP_total./(N_FP_total+N_TP_total);

save([foldername filesep files(1).name(1:end-6) '_STW_result.mat'],...
    'dt','dt_start','dt_stop','dt_start_array','dt_stop_array',...
    'dt_time','dt_time_start','dt_time_stop','dt_time_start_array','dt_time_stop_array',...
    'N_split','N_P','P','TPR','FDR','M','T','time_res','L');
%% plot and save figure
lw = 1;
fs = 8;
M_range = M(1,:);
T_range = T(:,1);
figure('Color',[1,1,1],'Units','inch','Position',[0,0,16,5]); hold on;
tiledlayout(2,5);

nexttile; hold on; colormap('vik');
imagesc(M(1,:),T(:,1),dt_start,'AlphaData',N_P>0);
%contour(M,T,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-5,5],'YLim',[0,T(end,1)],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('photons per time window, M');
ylabel('time window (\mus)');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_start(N_P == 0) = Inf;
[~,ix] = min(abs(dt_start),[],1,'linear');
plot(M(ix),T(ix),'.k');
%plot(M_range,M_range*mean(T(ix)./M(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(M(ix)./T(ix))*1E3));
c = colorbar;
c.Label.String = 'photon deviation of start';
title('photon deviation of start');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);
ax = gca;
text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','y','FontWeight','bold');

nexttile; hold on;
imagesc(M(1,:),T(:,1),dt_stop,'AlphaData',N_P>0);
%contour(M,T,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-5,5],'YLim',[0,T(end,1)],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('photons per time window, M');
ylabel('time window (\mus)');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_stop(N_P == 0) = Inf;
[~,ix] = min(abs(dt_stop),[],1,'linear');
plot(M(ix),T(ix),'.k');
%plot(M_range,M_range*mean(T(ix)./M(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(M(ix)./T(ix))*1E3));
c = colorbar;
c.Label.String = 'photon deviation of stop';
title('photon deviation of stop');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);

nexttile; hold on;
imagesc(M(1,:),T(:,1),dt,'AlphaData',N_P>0);
%contour(M,T,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,5],'YLim',[0,T(end,1)],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('photons per time window, M');
ylabel('time window (\mus)');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt(N_P == 0) = Inf;
[~,ix] = min(abs(dt),[],1,'linear');
plot(M(ix),T(ix),'.w');
%plot(M_range,M_range*mean(T(ix)./M(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(M(ix)./T(ix))*1E3));
c = colorbar;
c.Label.String = 'absolute photon deviation';
title('RMSD photon deviation');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);

nexttile; hold on;
imagesc(M(1,:),T(:,1),TPR,'AlphaData',~isnan(TPR) & N_P>0);
%contour(M,T,TPR,'LevelList',1,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,1],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'true positive rate';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('true positive rate');
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);

nexttile; hold on;
imagesc(M(1,:),T(:,1),N_P,'AlphaData',N_P>0);
contour(M,T,N_P,'LevelList',P,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,ceil(2.5*P/500)*500],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'number of bursts';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('number of detected bursts');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);
ax = gca;
text(ax.XLim(2)*0.525,ax.YLim(2)*0.95,sprintf('N_{true} = %i',round(P)),'Color','y','FontWeight','bold');

nexttile; hold on;
imagesc(M(1,:),T(:,1),dt_time_start*time_res*1E3,'AlphaData',N_P>0);
%contour(M,T,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-2,2],'YLim',[0,T(end,1)],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('photons per time window, M');
ylabel('time window (\mus)');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_time_start(N_P == 0) = Inf;
[~,ix] = min(abs(dt_time_start),[],1,'linear');
plot(M(ix),T(ix),'.k');
%plot(M_range,M_range*mean(T(ix)./M(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(M(ix)./T(ix))*1E3));
c = colorbar;
c.Label.String = 'time deviation of start (ms)';
title('time deviation of start (ms)');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(M(1,:),T(:,1),dt_time_stop*time_res*1E3,'AlphaData',N_P>0);
%contour(M,T,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-2,2],'YLim',[0,T(end,1)],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('photons per time window, M');
ylabel('time window (\mus)');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_time_stop(N_P == 0) = Inf;
[~,ix] = min(abs(dt_time_stop),[],1,'linear');
plot(M(ix),T(ix),'.k');
%plot(M_range,M_range*mean(T(ix)./M(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(M(ix)./T(ix))*1E3));
c = colorbar;
c.Label.String = 'time deviation of stop (ms)';
title('time deviation of stop (ms)');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(M(1,:),T(:,1),dt_time*time_res*1E3,'AlphaData',N_P>0);
%contour(M,T,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,2+9],'YLim',[0,T(end,1)],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('photons per time window, M');
ylabel('time window (\mus)');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_time(N_P == 0) = Inf;
[~,ix] = min(abs(dt_time),[],1,'linear');
plot(M(ix),T(ix),'.w');
%plot(M_range,M_range*mean(T(ix)./M(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(M(ix)./T(ix))*1E3));
c = colorbar;
c.Label.String = 'absolute time deviation (ms)';
title('RMSD time deviation (ms)');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(M(1,:),T(:,1),FDR,'AlphaData',~isnan(FDR) & N_P>0);
%contour(M,T,FDR,'LevelList',0,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,1],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'false discovery rate';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('false discovery rate');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(M(1,:),T(:,1),N_split,'AlphaData',N_P>0);
set(gca,'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'splitting frequency';
xlabel('photons per time window, M');
ylabel('time window (\mus)');
title('splitting frequency');
xlim([M_range(1),M_range(end)]); ylim([T_range(1),T_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

print(gcf,[foldername filesep files(1).name(1:end-6) '_STW.png'],'-dpng','-painters','-r300');
print(gcf,[foldername filesep files(1).name(1:end-6) '_STW.eps'],'-depsc','-painters');
print(gcf,[foldername filesep files(1).name(1:end-6) '_STW.pdf'],'-dpdf','-painters');
saveas(gcf,[foldername filesep files(1).name(1:end-6) '_STW.fig']);
delete(gcf);