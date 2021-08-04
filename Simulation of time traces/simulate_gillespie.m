function [MT,FileInfo,BurstPhotonNumbers] = simulate_gillespie(...
    Event_Rate,Duration,BurstDuration,I_event,I_background,time_res,Number_of_Levels,sigma_range)
% Input paramters:
% Event_Rate         - Event rate in Hz
% Duration           - Total duration in seconds
% BurstDuration      - The duration of a burst in milliseconds
% I_event            - The maximum burst intensity in kHz
% I_background       - The background signal in kHz
% time_res           - The time resolution in nanoseconds.
% Number_of_Levels   - Number of intensity levels (2N+1 steps in burst)
% sigma_range        - Plus/minus sigma to cover of Gaussian-shaped burst
%                      Set sigma_range = 0 to produce a square wave
%
% Output parameters:
% MT                 - The photon time stamps. (N_photonsx1) array
% FileInfo           - Stores information on the simulated dataset.
% BurstPhotonNumbers - The photon indices of photons detected during
%                      bursts. (N_burst x 1) cell array

% set default values for parameters
if ~exist('Event_Rate','var')
    Event_Rate = 5;
end
if ~exist('Duration','var')
    Duration = 1;
end
if ~exist('BurstDuration','var')
    BurstDuration = 1;
end
if ~exist('I_event','var')
    I_event = 20;
end
if ~exist('I_background','var')
    I_background = 1;
end
if ~exist('time_res','var')
    time_res = 12.5; % corresponding to 80 MHz
end
if ~exist('Number_of_Levels','var')
    Number_of_Levels = 3;
end
if ~exist('sigma_range','var')
    sigma_range = 1;
end

BurstDuration= BurstDuration/1000; %s
N_S=2*Number_of_Levels+1; %number of intensity states
time_res = time_res*1E-9;
I_event=I_event*1000; %Hz
I_background=I_background*1000; %Hz

I_out= I_event*exp(-linspace(-sigma_range,sigma_range,N_S).^2/2);
time_per_intensity = BurstDuration/N_S;
N_B = floor(Duration*Event_Rate); %number of bursts
ChangePoints=Event_Rate^(-1)*(0.5 + (0:N_B-1)); % evenly spaced timepoints of bursts in s
t_event = round(ChangePoints./time_res);

%%% Photon arrival time computation
% we simulate a fixed amount of time, based on a constant count rate

%disp('Simulating background.');
%tic
%draw the background photons
MT_Background = [];
if I_background > 0
    MT_Background = cumsum(ceil((exprnd(1/I_background,1,Duration*I_background)/time_res)))';
    while (MT_Background(end)<Duration/time_res) % we did not yet reach the end
        dt=ceil((exprnd(1/I_background))/time_res);
        MT_Background(end+1)=MT_Background(end)+dt;
    end
end
%toc

%draw event photons
%disp('Simulating events.');
%ll = fprintf('0 %%');
MT_Burst = zeros(1E6,1);
photon_index = 1;
% current_time  -   running variable to store the ablsoute time of the last
%                   photon
% burst_loop_time - running variable to store the relative time of the last
%                   photon with respect to the last intensity step
%tic
for i=1:numel(t_event)
    current_time=t_event(i)-round(BurstDuration/(2*time_res));
    for n=1:N_S
        burst_loop_time=0;
        while burst_loop_time<time_per_intensity/time_res
            dt=ceil((exprnd(1/(I_out(n))))/time_res);
            current_time=current_time+dt;
            burst_loop_time=burst_loop_time+dt;
            if current_time>t_event(i)-(BurstDuration/(2*time_res))...
                    +(((n)*time_per_intensity)/time_res) %if current time exceeds the intensity state time
                current_time=round(t_event(i)-(BurstDuration/(2*time_res))+(((n)*time_per_intensity)/time_res)); %reset the clock
                continue %go to the next loop without saving MT
            end
            % add photon
            if photon_index > numel(MT_Burst)
                MT_Burst = [MT_Burst;zeros(1E6,1)];
            end
            MT_Burst(photon_index) = current_time;
            photon_index = photon_index +1;
        end
    end
    % if mod(i,1000) == 0
    %   fprintf(repmat('\b',1,ll));
    %   ll = fprintf('%i %%\n',round(100*i/numel(t_event)));
    % end
end

%toc
MT_Burst(photon_index:end) = [];
MT = sort([MT_Background;MT_Burst]);%sort MTs of events and background in right time order
%if we are really unlucky a background photon and event photon are on the same timestamp,remove one photon
MT(find(diff(MT) == 0) +1 ) = [];

%remove photons exceeding the set duration
MT(MT>Duration/time_res)=[];

%disp('Finding ground truth burst starts and stops.');
%tic
% find photons in burst regions
BurstPhotonNumbers = cell(1,N_B);
for j = 1:N_B
   BurstPhotonNumbers{1,j} = find(abs(MT - t_event(j)) <= BurstDuration/(2*time_res));             
end
%toc

%%% Minimal set of meta data
FileInfo = struct;
FileInfo.ClockPeriod = time_res;
FileInfo.SyncPeriod = time_res; % the macrotime period
FileInfo.TACRange = time_res; % the microtime range
FileInfo.MeasurementTime = time_res*max(MT);
%%% store simulation input parameters
FileInfo.I_background= I_background; %Hz
FileInfo.I_event = I_event; %the highest intensity of an event
FileInfo.IntensityLevels = N_S; %number of intensity states
FileInfo.IntensityProfileSigma = sigma_range; % width of the normally-distributed intensity profile
FileInfo.Event_Rate = Event_Rate; % event rate in Hz
FileInfo.Duration = Duration; % duration in seconds
FileInfo.BurstDuration= BurstDuration; %s
FileInfo.NumberOfBursts = N_B; %number of bursts