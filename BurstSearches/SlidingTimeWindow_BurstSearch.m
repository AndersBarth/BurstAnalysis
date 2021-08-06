function [start,stop] = SlidingTimeWindow_BurstSearch(Photons,M,T,all_photons)
% M - Number of photons per time window
% T - time window length in microseconds
global FileInfo
if nargin < 4
    all_photons = true;
end

valid=(Photons(1+M-1:end)-Photons(1:end-M+1)) < T*1e-6/FileInfo.ClockPeriod;

% and find start and stop of bursts
if ~all_photons
    % legacy implementation before 07/2021
    % we take only the middle photon of a time window to belong to a burst
    start = find(valid(1:end-1)-valid(2:end)==-1) + 1 + floor(M/2); % +1 is necessary
    stop = find(valid(1:end-1)-valid(2:end)==1) + floor(M/2);
    %disp('Using legacy mode.');
else
    % instead, all photons of the time window are considered to be part of the burst
    start = find(valid(1:end-1)-valid(2:end)==-1) + 1; % +1 is necessary
    stop = find(valid(1:end-1)-valid(2:end)==1) + M - 1; % last photon and the M-1 following ones are included
    % make sure we have the same number of starts and stops
    if numel(start) < numel(stop)
        stop(1) = [];
    elseif numel(start) > numel(stop)
        start(end) = [];
    end
    % remove overlapping bursts which occur because "valid" only stores
    % the first photon of each burst
    while any(start(2:end)-stop(1:end-1) <= 0)
        overlap = find(start(2:end)-stop(1:end-1) <= 1); % next start is before current stop
        start(overlap+1) = [];
        stop(overlap) = [];
    end
end


if numel(start) < numel(stop)
    stop(1) = [];
elseif numel(start) > numel(stop)
    start(end) = [];
end