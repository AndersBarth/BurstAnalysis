clc
clear all
close all
format compact
Number_of_analysis=1;
%%% Initialisation and Intensity profiling
N_channels = 1;
I_background= 1000; %Hz
max_signal_to_noise_scale=100; %the highest intensity of an event
max_intensity= (max_signal_to_noise_scale)*I_background; %counts per event
I_event= I_background*max_signal_to_noise_scale; %currently I_event is kept constant at 20x snr
N_S=7; %number of intensity states
I_out= I_event*exp(-linspace(-1,1,N_S).^2/2);

Event_Rate = 5; % 5 events per second
Duration = 2000; % duration in seconds
burst_duration= 0.01; %s
time_per_intensity = burst_duration/N_S;
N_B=Duration*Event_Rate; %number of bursts
ChangePoints=Event_Rate^(-1)*(0.5 + (0:N_B-1)); % evenly spaced timepoints of bursts in s

time_res = 12.5*1E-9; % the time resolution of the system, here 12.5 nanosecond (80 MHz)
t_event= round(ChangePoints./time_res);
for iii=1:Number_of_analysis
    %% Photon arrival time computation
    % we simulate a fixed amount of time, based on a constant count rate
    current_time = 0; % last photon arrival time
    burst_loop_time=0; % to loop over the burst times
    %time_per_intensity=zeros(1,N_S);
    %burst_duration=zeros(1,N_S);
    %for i=1:numel(Spikes_sort)
    %   burst_duration(i)=(100)/1000; %in seconds
    %    time_per_intensity(i)=burst_duration(i)/(numel(x));
    %end
    %burst_duration(end+1)=0; %defined such that we can loop until the end of the simulation and not stop at last burst instantly
    MT = [];
    MT_Burst=[];
    MT_Background=[];
    PhotonNumbers=[];
    MT_save = cell(N_channels,1); % the "macrotime", i.e. the photon arrival time measured with respect to the start of the measurement
    
    disp('Simulating background.');
    tic
    %draw the background photons
    MT_Background = cumsum(ceil((exprnd(1/I_background,1,Duration*I_background)/time_res)));
    while (MT_Background(end)<Duration/time_res) % we did not yet reach the end
        dt=ceil((exprnd(1/I_background))/time_res);
        MT_Background(end+1)=MT_Background(end)+dt;
    end
    toc
    
    %draw event photons
    disp('Simulating events.');
    dtsave=[];
    tic
    for i=1:numel(t_event)
        current_time=t_event(i)-round(burst_duration/(2*time_res));
        for n=1:N_S
            burst_loop_time=0;
            while burst_loop_time<time_per_intensity/time_res
                dt=ceil((exprnd(1/(I_out(n))))/time_res);
                dtsave(end+1)=dt;
                if dt<1
                    dt=1;
                end
                current_time=current_time+dt;
                burst_loop_time=burst_loop_time+dt;
                if current_time>t_event(i)-(burst_duration/(2*time_res))+(((n)*time_per_intensity)/time_res) %if current time exceeds the intensity state time
                    current_time=round(t_event(i)-(burst_duration/(2*time_res))+(((n)*time_per_intensity)/time_res)); %reset the clock
                    continue %go to the next loop without saving MT
                end
                MT_Burst(end+1) = current_time;                
            end
        end
    end
    toc
    MT=sort([MT_Background,MT_Burst]);%sort MTs of events and background in right time order
    %if we are really unlucky a background photon and event photon are on the same timestamp,remove one photon
    MT(find(diff(MT) == 0) +1 ) = []; 
   
    %remove last photon such that every dataset ends in 1s
    MT(MT>Duration/time_res)=[];
    
    MT_save{1,1}=MT'; %Right format for PAM
    
    disp('Finding ground truth burst starts and stops.');
    tic
    % find photons in burst regions
    PhotonNumbers = cell(1,N_B);
    for j = 1:N_B
        PhotonNumbers{1,j} = find(abs(MT - t_event(j)) <= burst_duration/(2*time_res));             
    end
    toc
    % truncate to the shortest macrotime
    %MeasurementTime = min(cellfun(@max,MT_save));
    %for i = 1:numel(MT_save)
    %    MT_save{i} = MT_save{i}(MT_save{i} <= MeasurementTime);
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Below, additional data is defined %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % define "microtime" data (relates to the delay to the excitation pulse in
    % experiments with pulsed lasers)
    MI = cell(N_channels,1);
    MI_Bins = 100;
    for i = 1:N_channels
        MI{i,1} = randi(MI_Bins,size(MT_save{i,1}));
    end
    MI=MI';
    %%% Minimal set of meta data
    Info = struct;
    Info.ClockPeriod = time_res;
    Info.SyncPeriod = time_res; % the macrotime period
    Info.TACRange = time_res; % the microtime range
    Info.MI_Bins = MI_Bins; % number of microtime bins
    MeasurementTime =  max(cellfun(@max,MT_save));
    Info.MeasurementTime = time_res*MeasurementTime;
    % some other meta data that relates to imaging (but is required)
    Info.Pixels = 1;
    Info.Lines = 1;
    Info.LineTimes = [];
    %% data save for PAM
    %%% save the data
    MT=MT_save;
    out_path=''; %put your location you want to save on here
    FileName=sprintf('DataSet_%d.ppf',iii);
    full_path=fullfile(out_path,FileName);
    save(full_path,'MT','MI','Info');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if you want to save the photon data as well run the save command below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('PhotonBurstData.mat','PhotonNumbers')

%% plot
t = MT{1}; 
figure('Position',[100,100,1200,350]);
subplot(1,3,1); hold on;
plot(t*time_res);
xlabel('Photon number'); ylabel('Macrotime [s]');

subplot(1,3,2); hold on;
[h,time] = histcounts(t*time_res,0:1E-3:Duration);
plot(time(1:end-1),h./min(diff(time))/1000);
xlabel('Macrotime [s]'); ylabel('Count rate [kHz]');

subplot(1,3,3); hold on;
histogram(diff(t)*time_res*1E6);
xlabel('Interphoton time [µs]'); ylabel('Frequency');
set(gca,'YScale','log');