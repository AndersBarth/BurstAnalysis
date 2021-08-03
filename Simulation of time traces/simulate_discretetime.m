function [MT,FileInfo,BurstPhotonNumbers] = simulate_discretetime(...
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
N_B=floor(Duration*Event_Rate); %number of bursts
N_S=2*Number_of_Levels+1; %number of intensity states
TimePerState=(BurstDuration/(N_S)); %in seconds
time_res = time_res*1E-9;
Time_Event = Event_Rate^(-1)*(0.5 + (0:N_B-1)); %seconds
N = ceil(Duration/time_res); % number of timebins

I_event=I_event*1000; %Hz
I_background=I_background*1000; %Hz
%%% generate photons based on binomial distribution at each time interval
I_out_save = I_event*exp(-linspace(-sigma_range,sigma_range,7).^2/2);

Time_ChangePoints = zeros(N_B,N_S+1);
for i=1:N_B
    Time_ChangePoints(i,:)=Time_Event(i)-BurstDuration/2+(0:N_S)*TimePerState;
end

Time_ChangePoints_res = round(Time_ChangePoints/time_res); % round so we can use it as array indices

disp('Simulating photons...');
%tic
cr = zeros(N,1); % initialize array
for i = 1:size(Time_ChangePoints_res,1)
    for j = 1:N_S
        cr(Time_ChangePoints_res(i,j):Time_ChangePoints_res(i,j+1)) = I_out_save(j);
    end
end
cr = cr + I_background;%%% add background
p = time_res*cr; % probability to see a photon at any time point
MT = find(binornd(1,p)); % generate random numbers and find non-zero time points
%toc

disp('Finding ground truth burst starts and stops.');
%tic
% find photons in burst regions
BurstPhotonNumbers = cell(1,N_B);
Time_event_res = round(Time_Event./time_res);
for j = 1:N_B
    BurstPhotonNumbers{1,j} = find(abs(MT - Time_event_res(j)) <= BurstDuration/(2*time_res));             
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
