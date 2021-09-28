function [dsArray, dsFeatures, S] = directed_spectrum(data, areaList, fs, f,...
    dsVersion, opts)
% directed_spectrum
%    Estimates directed spectrum from multi-channel timeseries
%    data. This is an estimate of communication from one channel to another.
%
%    INPUTS
%    data: NxAxW array
%        Preprocessed (filtered, averaged, checked for saturation) data.
%        A is the # of areas. N = number of frequency points per
%        window. W = number of time windows.
%    areaList: cell array of strings
%        Ordered list of brain areas from which data was recorded. 
%    fs: scalar
%        Sampling rate (Hz) of the data
%    f: vector
%        Defines the (non-zero) frequencies at which the features should
%        be evaluted.
%    dsVersion: 2-element boolean vector
%        The first element indicates whether to calculate the 'full' model
%        for directed spectrum features; the second element indcates
%        whether to calculate the pairwise directed spectrum. Both
%        versions can be calculated together for comparison.
%    opts: structure
%    FIELDS
%      window: integer or vector
%          If an integer, gives length of Hamming subwindows used in
%          Welch's method for estimating a cross spectrum.
%          If a vector, indicates the window that should be used in
%          Welch's method. Can also be left as an empty array '[]' to use
%          the matlab default window size.
%      overlap: integer
%          Indicate the overlap between sucessive windows when using
%          Welch's method to estimate a cross spectrum. If left as empty
%          (i.e. []) then the matlab default is used.
%      parCores: integer
%          Indicates number of cores to use for parallel computing. If 0,
%          all tasks executed in serial. [Default: 0]
%    OUTPUTS
%    dsArray: PxFxW array
%        directed spectrum estimates. P iterates over
%        directed pairs of regions. F iterates over frequencies. W=number
%        of time windows.
%    dsFeatures: PxF array
%        Gives string labels describing the
%        features represented in dsArray. 
%    S: AxAxFxW complex array
%        Cross-power spectral density associated with each window. Where
%        A = # of areas; F = # of frequencies; W = # of windows.

% determine which frequencies to evaluate
f = [0, f];
F = numel(f);

[~,C,W] = size(data);
channelPairs = combnk(1:C, 2)';

% number of directed pairs
nCP = size(channelPairs, 2);

% generate list of strings for frequencies
fStrings = compose('%d', f)';

% initialize arrays to be filled in loop below
if dsVersion(1)
    dsArray12a = zeros(1, nCP, F, W, 'single');
    dsArray21a = zeros(1, nCP, F, W, 'single');
end
if dsVersion(2)
    dsArray12b = zeros(1, nCP, F, W, 'single');
    dsArray21b = zeros(1, nCP, F, W, 'single');
end
dsFeatures12 = cell(1, nCP, F);
dsFeatures21 = cell(1, nCP, F);
if nargout > 2
    S = zeros(C,C,F-1,W);
end
    
if opts.parCores, pp = parpool([2 opts.parCores]); end
a = tic;
for w = 1:W
    thisData = data(:, :, w);
    fprintf('Starting window %d: %2.1fs elapsed\n', w, toc(a))
    
    thisCpsd = cpsd(thisData, thisData, opts.window, opts.overlap, f, fs, 'mimo');
    thisCpsd = permute(thisCpsd, [2,3,1]);
    
    % Calculate full DS model
    if dsVersion(1)
        % decompose full CPSD
        [H, SIG] = wilson_sf(thisCpsd, fs);
        
        % calculate DS for all pairs
        parfor (cp = 1:nCP, opts.parCores)
            c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
            
            [ds12, ds21] = calc_ds(H, SIG, F, c1, c2);
            
            % fill in output array
            dsArray12a(1,cp,:,w) = ds12; dsArray21a(1,cp,:,w) = ds21;
        end
    end
    
    % Calculate pairwise model
    if dsVersion(2)
        parfor (cp = 1:nCP, opts.parCores)
            c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
            
            % decompose pairwise CPSD
            pairCpsd = thisCpsd([c1,c2],[c1,c2],:);
            [H, SIG] = wilson_sf(pairCpsd, fs);
            
            [ds12, ds21] = calc_ds(H, SIG, F, 1, 2);
            
            % fill in output array
            dsArray12b(1,cp,:,w) = ds12; dsArray21b(1,cp,:,w) = ds21;
        end
    end
    
    if nargout > 2
        % remove DC term from cpsd before saving
        S(:,:,:,w) = thisCpsd(:,:,2:end);
    end
    
    if w==1
        for cp = 1:nCP
            c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
            
            % save labels corresponing to each element of feature arrays above
            dsFeatures12(1, cp,:) = cellfun(@(x) [areaList{c1} '->' areaList{c2} ' ' x], ...
                fStrings, 'UniformOutput', false);
            dsFeatures21(1, cp,:) = cellfun(@(x) [areaList{c2} '->' areaList{c1} ' ' x], ...
                fStrings, 'UniformOutput', false);
        end
    end
    
end
if opts.parCores, delete(pp), end

%% combine values down here
if sum(dsVersion) == 2
    dsArray = {reshape([dsArray12a; dsArray21a], 2*nCP, F, W), ...
        reshape([dsArray12b; dsArray21b], 2*nCP, F, W)};
elseif dsVersion(1)
    dsArray = reshape([dsArray12a; dsArray21a], 2*nCP, F, W);
else
    dsArray = reshape([dsArray12b; dsArray21b], 2*nCP, F, W);
end

dsFeatures = reshape([dsFeatures12; dsFeatures21], 2*nCP, F);
end

function [dsxy, dsyx] = calc_ds(H, SIG, h, x, y)        
Hxy = H(x,y,:);
SIGyx = SIG(y,y)-(SIG(y,x)/SIG(x,x))*SIG(x,y); % partial covariance

Hyx = H(y,x,:);
SIGxy = SIG(x,x)-(SIG(x,y)/SIG(y,y))*SIG(y,x);

dsyx = zeros(1,h); dsxy = zeros(1,h);
for k = 1:h
    dsyx(k) = real(Hxy(:,:,k)*SIGyx*Hxy(:,:, k)');
    dsxy(k) = real(Hyx(:,:,k)*SIGxy*Hyx(:,:, k)');
end
end