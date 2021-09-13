%% run analysis on sliding time window
clear all;
%addpath(genpath('.'));
L = 5;
root_dir = 'M:\tnw\bn\cd\Shared\Anders\Simulations_BurstSearches';
files_SR = dir(root_dir);
for i = 3:numel(files_SR)
    files_BD = dir([root_dir filesep files_SR(i).name]);
        for ii = 3:numel(files_BD)
            files_BG = dir([root_dir filesep files_SR(i).name ...
                filesep files_BD(ii).name]);
            for iii = 3:numel(files_BG)
                files_SG = dir([root_dir filesep files_SR(i).name ...
                filesep files_BD(ii).name ...
                filesep files_BG(iii).name]);
                for iiii = 3:numel(files_SG)
                    current_folder = [root_dir filesep files_SR(i).name ...
                        filesep files_BD(ii).name ...
                        filesep files_BG(iii).name ...
                        filesep files_SG(iiii).name];
                    fprintf('Analyzing folder %s\n',current_folder);
                    InterphotonTime_parameter_scan(current_folder,L);
                    %CUSUM_parameter_scan(current_folder,L);
                    %SlidingTimeWindow_parameter_scan(current_folder,L);
                end
            end
        end
end
