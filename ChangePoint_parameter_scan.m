%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For the CUSUM, the interesting parameter seems to be the error rate %%%
%%% alpha. By default, it is set to 1/Nphotons, but the resulting low   %%%
%%% rate causes a high detection delay.                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChangePoint_parameter_scan(foldername,L)
% L - Minimum number of photons per burst
if nargin < 2
    L = 5;
end
time_res = 12.5E-9;% scan over threshold and false positive percentage
IT_range = 0.1:0.1:1;
alpha_range = logspace(-3,-0.05,25);
[IT,alpha] = meshgrid(IT_range,alpha_range);

% initialize averaged variables over all files
N_P_total = zeros(size(IT)); % number of detected bursts with at least L photons
dt_total = zeros(size(IT)); % mean photon deviation of start/stop
dt_start_total = zeros(size(IT)); % mean photon deviation of start
dt_stop_total = zeros(size(IT)); % mean photon deviation of stop
dt_start_array = cell(size(IT)); % array of all photon deviations of starts
dt_stop_array = cell(size(IT)); % array of all photon deviations of stops
dt_time_total = zeros(size(IT)); % mean time deviation of start/stop
dt_time_start_total = zeros(size(IT)); % mean time deviation of start
dt_time_stop_total = zeros(size(IT)); % mean time deviation of stop
dt_time_start_array = cell(size(IT)); % array of all time deviations of starts
dt_time_stop_array = cell(size(IT)); % array of all time deviations of stops
N_FP_total = zeros(size(IT)); % false positives
N_TP_total = zeros(size(IT)); % true positives
N_split_total = zeros(size(IT)); % splitting frequency
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
    dt = zeros(size(IT));
    dt_start = zeros(size(IT));
    dt_stop = zeros(size(IT));
    dt_time = zeros(size(IT));
    dt_start_time = zeros(size(IT));
    dt_stop_time = zeros(size(IT));
    
    N_FP = zeros(size(IT));
    N_P = zeros(size(IT));
    N_TP = zeros(size(IT));
    N_split = zeros(size(IT));
    for i = 1:size(IT,1)
        for j = 1:size(IT,2)
            [start,stop] = ChangePoint_BurstSearch(MT,...
                IT(i,j)*FileInfo.I_event/1000,alpha(i,j),false);
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
        %disp(i*size(IT,2)/numel(IT));
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

save([foldername filesep files(1).name(1:end-6) '_CUSUM_result.mat'],...
    'dt','dt_start','dt_stop','dt_start_array','dt_stop_array',...
    'dt_time','dt_time_start','dt_time_stop','dt_time_start_array','dt_time_stop_array',...
    'N_split','N_P','P','TPR','FDR','IT','alpha','time_res','L');
%% plot and save figure
lw = 1;
fs = 8;
IT_range = IT(1,:);
alpha = log10(alpha);
alpha_range = alpha(:,1);

figure('Color',[1,1,1],'Units','inch','Position',[0,0,16,5]); hold on;
tiledlayout(2,5);

nexttile; hold on; colormap('vik');
imagesc(IT(1,:),alpha(:,1),dt_start,'AlphaData',N_P>0);
%contour(IT,alpha,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-10,10],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_start(N_P == 0) = Inf;
[~,ix] = min(abs(dt_start),[],1,'linear');
plot(IT(ix),alpha(ix),'.k');
%plot(M_range,M_range*mean(alpha(ix)./IT(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(IT(ix)./alpha(ix))*1E3));
c = colorbar;
c.Label.String = 'photon deviation of start';
title('photon deviation of start');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);
ax = gca;
text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','y','FontWeight','bold');

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),dt_stop,'AlphaData',N_P>0);
%contour(IT,alpha,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-5,5],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_stop(N_P == 0) = Inf;
[~,ix] = min(abs(dt_stop),[],1,'linear');
plot(IT(ix),alpha(ix),'.k');
%plot(M_range,M_range*mean(alpha(ix)./IT(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(IT(ix)./alpha(ix))*1E3));
c = colorbar;
c.Label.String = 'photon deviation of stop';
title('photon deviation of stop');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),dt,'AlphaData',N_P>0);
%contour(IT,alpha,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,5],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt(N_P == 0) = Inf;
[~,ix] = min(abs(dt),[],1,'linear');
plot(IT(ix),alpha(ix),'.w');
%plot(M_range,M_range*mean(alpha(ix)./IT(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(IT(ix)./alpha(ix))*1E3));
c = colorbar;
c.Label.String = 'RMSD photon deviation';
title('RMSD photon deviation');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),TPR,'AlphaData',~isnan(TPR) & N_P>0);
%contour(IT,alpha,TPR,'LevelList',1,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,1],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'true positive rate';
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
title('true positive rate');
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),N_P,'AlphaData',N_P>0);
contour(IT,alpha,N_P,'LevelList',P,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,ceil(2.5*P/500)*500],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'number of bursts';
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
title('number of detected bursts');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);
ax = gca;
text(ax.XLim(2)*0.525,ax.YLim(2)*0.95,sprintf('N_{true} = %i',round(P)),'Color','y','FontWeight','bold');

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),dt_time_start*time_res*1E3,'AlphaData',N_P>0);
%contour(IT,alpha,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-2,2],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_time_start(N_P == 0) = Inf;
[~,ix] = min(abs(dt_time_start),[],1,'linear');
plot(IT(ix),alpha(ix),'.k');
%plot(M_range,M_range*mean(alpha(ix)./IT(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(IT(ix)./alpha(ix))*1E3));
c = colorbar;
c.Label.String = 'time deviation of start (ms)';
title('time deviation of start (ms)');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),dt_time_stop*time_res*1E3,'AlphaData',N_P>0);
%contour(IT,alpha,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[-2,2],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_time_stop(N_P == 0) = Inf;
[~,ix] = min(abs(dt_time_stop),[],1,'linear');
plot(IT(ix),alpha(ix),'.k');
%plot(M_range,M_range*mean(alpha(ix)./IT(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(IT(ix)./alpha(ix))*1E3));
c = colorbar;
c.Label.String = 'time deviation of stop (ms)';
title('time deviation of stop (ms)');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),dt_time*time_res*1E3,'AlphaData',N_P>0);
%contour(IT,alpha,N_P,'LevelList',1000,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,2+9],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
% plot guidelines at constant threshold
%plot(M_range,1E6*M_range/5000,'--');
%plot(M_range,1E6*M_range/10000,'--');
%plot(M_range,1E6*M_range/20000,'--');
dt_time(N_P == 0) = Inf;
[~,ix] = min(abs(dt_time),[],1,'linear');
plot(IT(ix),alpha(ix),'.w');
%plot(M_range,M_range*mean(alpha(ix)./IT(ix)),'--','LineWidth',lw);
%title(sprintf('Optimal threshold: %.2f kHz',mean(IT(ix)./alpha(ix))*1E3));
c = colorbar;
c.Label.String = 'RMSD time deviation (ms)';
title('RMSD time deviation (ms)');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),FDR,'AlphaData',~isnan(FDR) & N_P>0);
%contour(IT,alpha,FDR,'LevelList',0,'EdgeColor','y','LineWidth',lw);
set(gca,'CLim',[0,1],'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'false discovery rate';
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
title('false discovery rate');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

nexttile; hold on;
imagesc(IT(1,:),alpha(:,1),N_split,'AlphaData',N_P>0);
set(gca,'Box','on','LineWidth',lw,'FontSize',fs,'Layer','top');
c = colorbar;
c.Label.String = 'splitting frequency';
xlabel('rel. intensity threshold');
ylabel('input false positive rate, log \alpha');
title('splitting frequency');
xlim([IT_range(1),IT_range(end)]); ylim([alpha_range(1),alpha_range(end)]);
%ax = gca;
%text(ax.XLim(2)*0.775,ax.YLim(2)*0.95,sprintf('L = %i',L),'Color','w','FontWeight','bold');

print(gcf,[foldername filesep files(1).name(1:end-6) '_CUSUM.png'],'-dpng','-painters','-r300');
print(gcf,[foldername filesep files(1).name(1:end-6) '_CUSUM.eps'],'-depsc','-painters');
print(gcf,[foldername filesep files(1).name(1:end-6) '_CUSUM.pdf'],'-dpdf','-painters');
saveas(gcf,[foldername filesep files(1).name(1:end-6) '_CUSUM.fig']);
delete(gcf);