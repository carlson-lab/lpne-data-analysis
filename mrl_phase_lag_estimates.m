function [lagXY, mrl] = mrl_phase_lag_estimates(x, y, f, fs, filtRatio, lags)
% Estimate lag via mean resultant length of phase offsets
% Estimates the time lag between signals x and y associated with specified
% frequency bands. Estimates are given for all frequency bands given by f.
% Negative values indicate signal x leads signal y. Lag is calculated by:
%    1.calculating instantaneous phase of each signal via Hilbert transform
%    2. For a single time lag calculate instantaneous phase offset
%    (difference) between x and y
%    3. evaluate mean resultant length of all phase offset values in series 
%    4. repeat steps 2/3 for all time lags of interest
%    5. choose lag with greates mean resultant length.
% DO NOT USE WITH WINDOWS SHORTER THAN 1s!
%  
% INPUTS 
% x, y: Signals. Each should be a NxW matrix of time windows (N: Time
%    points per window; W: number of windows).
% f: List of frequency band boundaries. 2xF matrix (F: Number of bands)
% fs: Sampling rate (Hz)
% filtRatio: (Optional) Fraction of the sampling frequency to use
%     as the filter order for frequency band filter. Can also be
%     thought of as the fraction of a 1s window that is lost due to
%     filtering delay.
% lags: (Optional) Vector of time lags to consider (ms) Default is
%     -150:150 ms
%
% OUTPUTS
% lagXY: FxW matrix of estimated lags where F is frequency and W is
%     window length
if nargin < 6
    lags = -100:100; % lag in ms
    if nargin < 5
        filtRatio = 1/4;
    end
end
lagPts = round(lags * fs / 1e3); % lag in number of timepoints

if size(f, 1) ~= 2
    if size(f, 1) == 1
        % convert list of frequencies to bands by assuming 1Hz bandwidths
        f = bsxfun(@plus, f, [-0.5; 0.5]);
    elseif size(f, 2) == 1
        f = bsxfun(@plus, f.', [-0.5; 0.5]);
    else
        error('Dimensions of f are incorrect!\n')
    end
end

W = size(x,2); % number of windows
F = size(f,2); % number of frequencies
               % initialize matrix of lag times
lagXYPts = zeros(F,W);
mrl = zeros(F,numel(lags),W);

t = tic;
for n = 1:F
    
    % Filter signals
    b = fir1(fs*filtRatio, f(:,n) * 2/fs);
    delay = ceil(mean(grpdelay(b,1)));
    xFiltered = filter(b, 1, double(x));
    xFiltered(1:delay,:) = [];
    yFiltered = filter(b, 1, double(y));
    yFiltered(1:delay,:) = [];

    % find estimated lag for each frequency
    [lagXYPts(n, :), mrl(n, :, :)] = lag_estimate(xFiltered, yFiltered, lagPts);

    fprintf('MRL calculations for %dHz to %dHz complete: %.1fs elapsed\n', ...
        f(1,n), f(2,n), toc(t));
end

% convert to ms for lag time
lagXY = lagXYPts * 1e3 / fs;
end

function [lagEst, allMRLs] = lag_estimate(x, y, lags)
% Calculates lag estimate via mean resultant length of instantaneous phase
% offset.
%
% INPUTS
% x,y: Filtered signals. Each should be a NxW matrix of time windows (N: Time
%        points per window; W: number of windows).
% lags: lag in number of timepoints
% 
% OUPUTS
% lagEst: estimated lag based on maximum mean resultant of offset

L = numel(lags);
W = size(x, 2);

% get instantaneous phase
xPhase = angle(hilbert(x));
yPhase = angle(hilbert(y));

% calculate mean resultant length for each lag
meanResultant = zeros(L, W);

% for each lag, shift X and Y by lag in each direction so they are
% matched. Take offset of x-y
parfor l = 1:L
    thisLag = lags(l);
    if thisLag < 0
        thisLag = -thisLag;
        % correct for x lags y
        xPhaseShifted = xPhase(thisLag+1 : end, :);
        yPhaseShifted = yPhase(1 : end - thisLag, :);
    else
        % correct for x lead y
        xPhaseShifted = xPhase(1 : end - thisLag, :);
        yPhaseShifted = yPhase(thisLag+1 : end, :);
    end
    
    offset = xPhaseShifted - yPhaseShifted;
    meanResultant(l,:) = mean( exp(1i .* offset));
end
allMRLs = abs(meanResultant);

% pick the lag that maximizes mean of the resultant
[~, estIdx] = max(allMRLs);

lagEst = lags(estIdx);
end

