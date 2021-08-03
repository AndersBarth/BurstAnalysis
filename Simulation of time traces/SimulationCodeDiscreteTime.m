iii = 1;
Event_Rate = 5; % 5 events per second
Duration = 60; % duration in seconds
BurstDuration= 0.010; %s
N_B=floor(Duration*Event_Rate); %number of bursts
N_S=7; %number of intensity states
TimePerState=(BurstDuration/(N_S)); %in seconds
sigma_range = 1; % how much sigma to cover of the normally-distributed intensity profile
Time_event = Event_Rate^(-1)*(0.5 + (0:N_B-1)); %seconds

res = 12.5E-9; % time resolution in seconds
N = ceil(Duration/res); % number of timebins

I_event=20000; %Hz
I_background=1000; %Hz
%%% generate photons based on binomial distribution at each time interval
I_out_save = zeros(5,N_S);
for i= 1:N_B
    I_out_save(i,:)= I_event*exp(-linspace(-sigma_range,sigma_range,7).^2/2);
end

Time_ChangePoints = zeros(N_B,N_S+1);
for i=1:N_B
    Time_ChangePoints(i,:)=Time_event(i)-BurstDuration/2+(0:N_S)*TimePerState;
end

Time_ChangePoints_res = round(Time_ChangePoints/res); % round so we can use it as array indices

chunksize = 1E8;
N_chunks = ceil(N./chunksize);
t = cell(N_chunks,1);
for c = 1:N_chunks
    if c == N_chunks % last bin, go up to end
        from_to = [(c-1)*chunksize+1,N];
    else
        from_to = [(c-1)*chunksize+1,c*chunksize];
    end
    cr = zeros(from_to(2)-from_to(1)+1,1); % initialize array
    valid = Time_ChangePoints_res(:,1) >= from_to(1) & Time_ChangePoints_res(:,1) <= from_to(2);
    cp_temp = Time_ChangePoints_res(valid,:) - from_to(1);
    for i = 1:size(cp_temp,1)
        for j = 1:N_S
            cr(cp_temp(i,j):cp_temp(i,j+1)) = I_out_save(i,j);
        end
    end
    cr = cr + I_background;%%% add background
    p = res*cr; % probability to see a photon at any time point
    t{c} = find(binornd(1,p)) + from_to(1); % generate random numbers and find non-zero time points
    disp(c/N_chunks);
end
t = vertcat(t{:});

disp('Finding ground truth burst starts and stops.');
tic
% find photons in burst regions
PhotonNumbers = cell(1,N_B);
Time_event_res = round(Time_event./res);
for j = 1:N_B
    PhotonNumbers{1,j} = find(abs(t - Time_event_res(j)) <= BurstDuration/(2*res));             
end
toc

MT={t};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Below, additional data is defined %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define "microtime" data (relates to the delay to the excitation pulse in
% experiments with pulsed lasers)
MI = cell(1,1);
MI_Bins = 100;
for i = 1:1
    MI{i,1} = randi(MI_Bins,size(t));
end
MI=MI';
%%% Minimal set of meta data
Info = struct;
Info.ClockPeriod = res;
Info.SyncPeriod = res; % the macrotime period
Info.TACRange = res; % the microtime range
Info.MI_Bins = MI_Bins; % number of microtime bins
Info.MeasurementTime = Duration;
% some other meta data that relates to imaging (but is required)
Info.Pixels = 1;
Info.Lines = 1;
Info.LineTimes = [];
%%% store simulation input parameters
FileInfo.I_background= I_background; %Hz
FileInfo.I_event = I_event; %the highest intensity of an event
FileInfo.IntensityLevels = N_S; %number of intensity states
FileInfo.IntensityProfileSigma = sigma_range; % width of the normally-distributed intensity profile
FileInfo.Event_Rate = Event_Rate; % event rate in Hz
FileInfo.Duration = Duration; % duration in seconds
FileInfo.BurstDuration= BurstDuration; %s
FileInfo.NumberOfBursts = N_B; %number of bursts
%% data save for PAM
%%% save the data
out_path=''; %put your location you want to save on here
FileName=sprintf('DataSet_%d.ppf',iii);
full_path=fullfile(out_path,FileName);
save(full_path,'MT','MI','FileInfo','PhotonNumbers');
%% plot
if 0
    dt = diff(t); % calculate interphoton times
    figure('Position',[100,100,1200,350]);
    subplot(1,3,1); hold on;
    plot(t*res);
    xlabel('Photon number'); ylabel('Macrotime [s]');

    subplot(1,3,2); hold on;
    [h,time] = histcounts(t*res,0:1E-3:Duration);
    plot(time(1:end-1),h./min(diff(time))/1000);
    xlabel('Macrotime [s]'); ylabel('Count rate [kHz]');

    subplot(1,3,3); hold on;
    histogram(dt*res*1E6);
    xlabel('Interphoton time [µs]'); ylabel('Frequency');
    set(gca,'YScale','log');
end
%% plot
if 0
    figure()
    plot((1:N)*res,cr/1000,'LineWidth',1.8)
    hold on
    bins=0:0.0133333333333333333333333333333333333333333:1;
    counts=histcounts(MT{1}*res,bins);
    normalcounts=(counts/max(counts))*21;
    x=bins(1:(end)-1);
    stairs(x,normalcounts,':','LineWidth',2.5)
    ylim([0 26])

    xlabel('Time in seconds','FontSize',15)
    ylabel('Intensity in kHz','FontSize',15)
    legend('Intended intensity profile' ,'Simulated Photon arrival times intensity profile','FontSize',11)
end