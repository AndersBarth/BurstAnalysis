%% run analysis on sliding time window
clear all;
L = 5;
files_SR = dir('Simulations');
for i = 3:numel(files_SR)
    files_BD = dir(['Simulations' filesep files_SR(i).name]);
        for ii = 3:numel(files_BD)
            files_BG = dir(['Simulations' filesep files_SR(i).name ...
                filesep files_BD(ii).name]);
            for iii = 3:numel(files_BG)
                files_SG = dir(['Simulations' filesep files_SR(i).name ...
                filesep files_BD(ii).name ...
                filesep files_BG(iii).name]);
                for iiii = 3:numel(files_SG)
                    current_folder = ['Simulations' filesep files_SR(i).name ...
                        filesep files_BD(ii).name ...
                        filesep files_BG(iii).name ...
                        filesep files_SG(iiii).name];
                    fprintf('Analyzing folder %s\n',current_folder);
                    SlidingTimeWindow_parameter_scan(current_folder,L);
                end
            end
        end
end
