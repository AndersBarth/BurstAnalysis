%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function for CUSUM based burst search                              %%%%
%%% Based on: Zhang, K. & Yang, H. Photon-by-photon determination of   %%%%
%%% emission bursts from diffusing single chromophores.                %%%%
%%% J Phys Chem B 109, 21930-21937 (2005).                             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [START,STOP] = CUSUM_BurstSearch(Photons,IB,IT)
global FileInfo
START = [];
STOP = [];
% convert photons to delay times
dt = diff(Photons);
if any(mod(dt,1)~=0)
    disp('Warning: Non-integer macrotimes found. Interphoton times will be rounded.');
    dt = round(dt);
end
% parameters (all count rates are given in kHz and need to be converted to
% the used clock period)
IB = IB*1E3*FileInfo.ClockPeriod; % background count rate per clock period
IT = IT*1E3*FileInfo.ClockPeriod; % threshold intensity
% alternative parameterization based on molecular brightness to determine
% the threshold
% I0 = I0*1E3*FileInfo.ClockPeriod; % intensity of molecule at the center of PSF (molecular brightness)
% IT = I0*exp(-2)+IB; % threshold intensity at 1/e^2

% error rates
alpha = 1/numel(Photons);
beta = 0.05;
% calculate the expectation value of the log likelihood ratio
x = 0:1:max([max(dt), ceil(-log(1E-3)/IB)]);
fB = exp(-x*IB); fB = fB./sum(fB); fB(fB==0) = eps;
fT = exp(-x*IT); fT = fT./sum(fT); fT(fT==0) = eps;
lambda = log(fT)-log(fB);
%mlambdaB = sum(fB.*lambda);
mlambdaT = sum(fT.*lambda);
% define CUSUM threshold h
h = -log(alpha*log(alpha^(-1))/(3*(mlambdaT+1)^2)); % eq. 4
% define SPRT thresholds A and B
B = beta/(1-alpha);
A = (1-beta)/alpha;

start_next = 1; % of the next burst
stop = 1; % stop of the previous burst, start searching from here
while start_next < numel(dt) % we have not reached the end
    % find the first edge using CUSUM
    start = CUSUM(dt,h,fB,fT,stop);
    
    % estimate the end of burst using SPRT
    stop_est = SPRT(dt,A,B,fB,fT,start);
    if stop_est > start+5 % require an offset of at least 5 photons
        % find the next edge using CUSUM
        start_next = CUSUM(dt,h,fB,fT,stop_est);

        % do backwards CUSUM to refine the end of the previous burst
        stop = bCUSUM(dt,h,fB,fT,start,start_next,10);

        if stop < stop_est && stop > start
            if isempty(START) || START(end) ~= start
                START(end+1,1) = start;
                STOP(end+1,1) = stop;
            else % sometimes, the algorithm gets stuck
                % move on
                stop = start + 10;
                start_next = stop; % to trigger exit condition
            end
        end
    else % sometimes, the algorithm gets stuck
        % move on
        stop = start + 10;
        start_next = stop; % to trigger exit condition
    end
end

function ix = CUSUM(dt,h,fB,fT,ix_start)
% find the first edge using CUSUM
if nargin < 5
    ix_start = 1;
end
ix = ix_start-1; % go one back because we increase the counter in the while loop before evaluating the S function
S = 0;
while S < h && ix < numel(dt)
    ix = ix + 1;
    S = max([S+log(fT(dt(ix)+1))-log(fB(dt(ix)+1)),0]);
end

function ix = SPRT(dt,A,B,fB,fT,ix_start)
% estimate the burst end using SPRT
if nargin < 6
    ix_start = 1;
end
ix = ix_start;
LAMBDA = fT(dt(ix_start)+1)/fB(dt(ix_start)+1);
while LAMBDA > B && ix < numel(dt)
    ix = ix+1;
    LAMBDA = LAMBDA*fT(dt(ix)+1)/fB(dt(ix)+1);
    if LAMBDA >= A
        LAMBDA = A;
    elseif LAMBDA <= B
        LAMBDA = B;
    end
end

function ix = bCUSUM(dt,h,fB,fT,ix_start,ix_next,offset)
if nargin < 7
    offset = 10; % offset necessary because if we start in the burst, the threshold is crossed in the beginning already.
end
dt_b = dt(ix_next-offset:-1:ix_start);
ix = ix_next-(offset-1)-CUSUM(dt_b,h,fB,fT)+1; % plus one needed to correctly invert the index
