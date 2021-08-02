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
N_B=floor(Duration*Event_Rate); %number of bursts
N_S=2*Number_of_Levels+1; %number of intensity states
TimePerState=(BurstDuration/(N_S)); %in seconds
time_res = time_res*1E-9;
Time_Event = Event_Rate^(-1)*(0.5 + (0:N_B-1)); %seconds
