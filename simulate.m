%%% Perform simulation over a wide range of conditions using the
%%% Gillespie algorithm implementation
addpath(genpath('.'));
%%% Basic settings
%   Simulate a total of 10000 bursts at 1 burst/second
Duration = 1000; % simulation duration in s
repeats = 10; % number of simulations per condition
Event_Rate = 1; % Hz
Time_Res = 12.5; % ns, default time resolution of 80 MHz

%%% Burst shape
Burst_Duration = [0.1,1,10,100]; % in ms
Number_of_Levels = 3; % 3 intensity levels for all simulations
sigma_range = [0,1,2]; % sample square wave and over range from -2/-1 to +1/+2 sigma


%%% Signal to noise ratio
I_event = [10,20,50,100]; % signal intensities to sample in kHz
I_background = [0,1,5,10]; % background count rates to sample

%%% Run over all pairs of parameters
% rmdir('Simulations');
%mkdir('Simulations');
N_sim = numel(sigma_range)*numel(Burst_Duration)*numel(I_background)*numel(I_event);
fprintf('Starting simulation...\n');
for i = 1:numel(sigma_range)
    for ii = 1:numel(Burst_Duration)
        for iii = 1:numel(I_background)
            for iiii = 1:numel(I_event)
                current_folder = get_folder_for_simulation(sigma_range(i),...
                    Burst_Duration(ii),I_background(iii),I_event(iiii));
                mkdir(current_folder);
                for r = 1:repeats
                    %%% do simulation
                    [MT,FileInfo,BurstPhotonNumbers] = simulate_gillespie(...
                        Event_Rate,Duration,Burst_Duration(ii),I_event(iiii),...
                        I_background(iii),Time_Res,Number_of_Levels,...
                        sigma_range(i));
                    %%% generate filename
                    filename = sprintf('Sim_SigR%i_Dur%g_BG%i_I%i_%i.ppf',...
                        sigma_range(i),Burst_Duration(ii),I_background(iii),...
                        I_event(iiii),r);
                    %%% save                
                    save([current_folder filesep filename],...
                        'MT','FileInfo','BurstPhotonNumbers');
                end
                fprintf('Finished simulation of condition %i of %i\n',(i-1)*numel(Burst_Duration)*numel(I_event)*numel(I_background)+ ...
                    (ii-1)*numel(I_event)*numel(I_background) + ...
                    (iii-1)*numel(I_background) + ...
                    iiii,N_sim);   
            end
        end
    end
end



function folder = get_folder_for_simulation(a,b,c,d)
folder = ['Simulations' filesep sprintf('SigmaRange_%i',a) ...
            filesep sprintf('BurstDuration_%g',b) ...
            filesep sprintf('Background_%ikHz',c) ...
            filesep sprintf('Signal_%ikHz',d)];
end