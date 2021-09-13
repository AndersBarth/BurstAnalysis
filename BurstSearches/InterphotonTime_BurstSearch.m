function [start,stop] = InterphotonTime_BurstSearch(Photons,m,T)
global FileInfo
% Seidel-Type burstsearch based on interphoton time and Lee Filter
% m - the smoothing window size (2m+1 photons)
% T - the interphoton time threshold in microseconds
if m>1
    % Smooth the interphoton time trace
    dT =[Photons(1);diff(Photons)];
    dT_m = zeros(size(dT,1),size(dT,2));
    dT_s = zeros(size(dT,1),size(dT,2));
    % Apply Lee Filter with window 2m+1
    sig_0 = std(dT); %%% constant filter parameter is the noise variance, set to standard devitation of interphoton time
    dT_cumsum = cumsum(dT);
    dT_cumsum = [0; dT_cumsum];
    for i = m+1:numel(dT)-m
        dT_m(i) = (2*m+1)^(-1)*(dT_cumsum(i+m+1)-dT_cumsum(i-m));
    end
    dT_sq_cumsum = cumsum((dT-dT_m).^2);
    dT_sq_cumsum = [0;dT_sq_cumsum];
    for i = 2*m:numel(dT)-2*m
        dT_s(i) = (2*m+1)^(-1)*(dT_sq_cumsum(i+m+1)-dT_sq_cumsum(i-m));
    end

    %filtered data
    dT_f = dT_m + (dT-dT_m).*dT_s./(dT_s+sig_0.^2);

    % threshold
    valid = dT_f < T & dT_f > 0;

elseif m == 1
    % threshold
    valid = [Photons(1); diff(Photons)] < T;
end
% and find start and stop of bursts
start = find(valid(1:end-1)-valid(2:end)==-1);
stop = find(valid(1:end-1)-valid(2:end)==1);

if numel(start) < numel(stop)
    stop(1) = [];
elseif numel(start) > numel(stop)
    start(end) = [];
end