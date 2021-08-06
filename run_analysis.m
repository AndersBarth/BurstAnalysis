%% run analysis on sliding time window
files = dir('Simulations/SigmaRange_0/BurstDuration_1/Background_1kHz/Signal_20kHz');

for i = 3:numel(files)
    SlidingTimeWindow_parameter_scan(files(i).name,L);
end