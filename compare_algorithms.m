clear all
dur = 2000;
event_rate = 5;
burst_dur = 100;

[MT_g,I_g,P_g] = simulate_gillespie(event_rate,dur,burst_dur);
N_phot = 0;
for i = 1:dur
    [MT_d{i},I_d{i},P_d{i}] = simulate_discretetime(event_rate,1,burst_dur);
    MT_d{i} = MT_d{i} + ceil((i-1)/12.5E-9);
    P_d{i} = cellfun(@(x) x + N_phot,P_d{i},'UniformOutput',false);
    N_phot = N_phot + numel(MT_d{i});
    i
end
MT_d = vertcat(MT_d{:});
P_d = horzcat(P_d{:});
%%
figure;hold on;
plot(histcounts(MT_d*12.5E-9,0:1E-3:dur));
plot(histcounts(MT_g*12.5E-9,0:1E-3:dur));

%% compare starts
start_gt_d = MT_d(cellfun(@(x) x(1), P_d));
stop_gt_d = MT_d(cellfun(@(x) x(end), P_d));
start_gt_g = MT_g(cellfun(@(x) x(1), P_g));
stop_gt_g = MT_g(cellfun(@(x) x(end), P_g));

% true start stops at
event_times = event_rate^(-1)*(0.5 + (0:floor(event_rate*dur)-1)'); 
start_true = event_times - burst_dur/1000/2;
stop_true = event_times + burst_dur/1000/2;
%% plotting
lw = 0.75;
fs = 8;
hf = figure('Units','inch','Position',[0,0,8,4]);
time_res = 12.5E-9*1E6; % mus
ax(1) = subplot(2,3,1);
hold on;
histogram(time_res*start_gt_g-start_true*1E6,time_res*(0:1E3:5E4),'EdgeColor','none');
histogram(time_res*start_gt_d-start_true*1E6,time_res*(0:1E3:5E4),'EdgeColor','none');
xlabel('time [µs]');
ylabel('#');
title('deviation of start from GT');
legend({'Gillespie','discrete time'},'EdgeColor','none','Color','none');

ax(2) = subplot(2,3,2);
hold on;
histogram(-time_res*stop_gt_g+stop_true*1E6,time_res*(0:1E3:5E4),'EdgeColor','none');
histogram(-time_res*stop_gt_d+stop_true*1E6,time_res*(0:1E3:5E4),'EdgeColor','none');
xlabel('time [µs]');
ylabel('#');
title('deviation of stop from GT');

ax(3) = subplot(2,3,3);
hold on;
bins = -200.5:5:200.5;
histogram(time_res*(stop_gt_g+start_gt_g)/2-event_times*1E6,bins,'EdgeColor','none');
histogram(time_res*(stop_gt_d+start_gt_d)/2-event_times*1E6,bins,'EdgeColor','none');
xlabel('time [µs]');
ylabel('#');
title('deviation of mean from GT');

ax(4) = subplot(2,3,4);
hold on;
bins = -500.5:10:500.5;
histogram(time_res*(start_gt_g-start_gt_d),bins,'EdgeColor','none');
xlabel('time [µs]');
ylabel('#');
title('difference between starts');

ax(5) = subplot(2,3,5);
hold on;
histogram(time_res*(stop_gt_g-stop_gt_d),bins,'EdgeColor','none');
xlabel('time [µs]');
ylabel('#');
title('difference between stops');

ax(6) = subplot(2,3,6);
hold on;
bins = linspace(0.99,1,100);
histogram((stop_gt_g-start_gt_g)*time_res/(1000*burst_dur),bins,'EdgeColor','none');
histogram((stop_gt_d-start_gt_d)*time_res/(1000*burst_dur),bins,'EdgeColor','none');
xlabel('relative burst length');
ylabel('#');
title('burst length');

%set(ax(1:2),'YScale','log');
ax(1).XLim(1) = 0;
ax(2).XLim(1) = 0;
set(ax,'FontSize',fs,'LineWidth',lw,'Box','on','FontName','Arial');

print(gcf,'comparison_of_algorithms.png','-dpng','-painters','-r300');
print(gcf,'-dpdf','comparison_of_algorithms.pdf','-painters');