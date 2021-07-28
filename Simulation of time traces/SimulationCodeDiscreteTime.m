Number_of_Traces = 1;
Event_Rate = 5; % 5 events per second
Duration = 1; % duration in seconds
BurstDuration= 0.010; %s
N_B=Duration*Event_Rate; %number of bursts
N_S=7; %number of intensity states
TimePerState=(BurstDuration/(N_S)); %in seconds

Time_event = Event_Rate^(-1)*(0.5 + (0:N_B-1)); %seconds

res = 100E-9; % time resolution in seconds
N = ceil(Duration/res); % number of timebins

I_event=20000; %Hz
I_background=1000; %Hz
%%% generate photons based on binomial distribution at each time interval
for iii=1:Number_of_Traces  
    I_out_save = zeros(5,N_S);
    for i= 1:N_B
        I_out_save(i,:)= I_event*exp(-linspace(-1,1,7).^2/2);
    end
    
    for i=1:N_B
        Time_ChangePoints(i,:)=Time_event(i)-BurstDuration/2+(0:N_S)*TimePerState;
    end
    
    Time_ChangePoints_res=Time_ChangePoints/res;
    cr = zeros(N,1); % initialize array
    Time_ChangePoints_res = round(Time_ChangePoints_res); % round so we can use it as array indices
    for i = 1:N_B
        for j = 1:N_S
            cr(Time_ChangePoints_res(i,j):(Time_ChangePoints_res(i,j+1))) = I_out_save(i,j);
        end
    end
    
    cr = cr + I_background;%%% add background
    p = res*cr; % probability to see a photon at any time point
    t = find(binornd(1,p,N,1)); % generate random numbers and find non-zero time points
    dt = diff(t); % calculate interphoton times
    MT={t};
    
    if exist('PhotonNumbersSecondApproach100_0ms','var')==0
        PhotonNumbersSecondApproach100_0ms=cell(1000,N_B);
    end
    for i=1:numel(t)
        for ii=1:N_B
            for jj=1:N_S
                if t(i)>=Time_ChangePoints_res(ii,jj) && t(i)<=Time_ChangePoints_res(ii,jj+1)
                    %PhotonNumbersSecondApproach100_0ms{iii,ii}(end+1)=i;
                end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Below, additional data is defined %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % truncate to the shortest macrotime
    MeasurementTime = min(cellfun(@max,MT));
    for i = 1:numel(MT)
        MT{i} = MT{i}(MT{i} <= MeasurementTime);
    end
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
    Info.MeasurementTime = res*MeasurementTime;
    % some other meta data that relates to imaging (but is required)
    Info.Pixels = 1;
    Info.Lines = 1;
    Info.LineTimes = [];
    %% data save for PAM
    %%% save the data
    %out_path='C:\Users\odink\Documents\Uni\BEP\Matlab\Simulations\Intensity 20x, burst vary\Approach2\Sliding Time Window\Intensity_20x_background__burstduration_100.0ms';
    %[FileName, PathName] = uiputfile({'*.ppf','PAM photon file'},'Select file name','data.ppf');
    %FileName=sprintf('DataSet_%d.ppf',iii);
    %full_path=fullfile(out_path,FileName);
    %save(full_path,'MT','MI','Info');
end
%save('100.0msBurstsPhotonDataSecondApproach.mat','PhotonNumbersSecondApproach100_0ms');
%clear all
%% plot
figure('Position',[100,100,1200,350]);
subplot(1,3,1); hold on;
plot(t*res);
xlabel('Photon number'); ylabel('Macrotime [s]');

subplot(1,3,2); hold on;
[h,time] = histcounts(t*res,100);
plot(time(1:end-1),h./min(diff(time))/1000);
xlabel('Macrotime [s]'); ylabel('Count rate [kHz]');

subplot(1,3,3); hold on;
histogram(dt*res*1E6);
xlabel('Interphoton time [µs]'); ylabel('Frequency');
set(gca,'YScale','log');
%% plot
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
