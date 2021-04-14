function [ldArray, ldFeatures, S] = linear_directionality(data, areaList, fs, f,...
    ldVersion, opts)
% linear_directionality
%    Estimates linear directionality spectrum from multi-channel timeseries
%    data. This is an estimate of communication from one channel to another.
%
%    INPUTS
%    data: NxAxW array
%        Preprocessed (filtered, averaged, checked for saturation) data.  A is
%        the # of areas. N=number of frequency points per
%        window. W=number of time windows.
%    areaList: cell array of strings
%        Ordered list of brain areas from which data was recorded   . 
%    fs: scalar
%        Sampling rate (Hz) of the data
%    f: vector
%        Defines the (non-zero) frequencies at which the features should be
%        evaluted.
%    ldVersion: 2-element boolean vector
%        The first element indicates whether to calculate the 'full' model for
%        linear directionality features; the second element indcates
%        whether to calculate them pairwise. Both versions can be
%        calculated together for comparison.
%    opts: structure
%    FIELDS
%      window: integer or vector
%          If an integer, gives length of Hamming subwindows used in
%          Welch's method for estimating a power spectrum.
%          If a vector, indicates the window that should be used in Welch's
%          method. Can also be left as an empty array '[]' to use the
%          matlab defaule window size.
%      overlap: integer
%          Indicate the overlap between sucessive windows when using
%          Welch's method to estimate a power spectrum. If left as empty
%          (i.e. []) then the matlab default is used.
%      parCores: integer
%          Indicates number of cores to use for parallel computing. If 0,
%          all tasks executed in serial. [Default: 0]
%    OUTPUTS
%    ldArray: PxFxW array
%        linear directionality estimates. P iterates over
%        directed pairs of regions. F iterates over frequencies. W=number
%        of time windows.
%    ldFeatures: PxF array
%        Gives string labels describing the
%        features represented in ldArray. 
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
if ldVersion(1)
    ldArray12a = zeros(1, nCP, F, W, 'single');
    ldArray21a = zeros(1, nCP, F, W, 'single');
end
if ldVersion(2)
    ldArray12b = zeros(1, nCP, F, W, 'single');
    ldArray21b = zeros(1, nCP, F, W, 'single');
end
ldFeatures12 = cell(1, nCP, F);
ldFeatures21 = cell(1, nCP, F);
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
    
    % Calculate full directionality model
    if ldVersion(1)
        % decompose full CPSD
        [H, SIG] = wilson_sf(thisCpsd, fs);
        
        % calculate LSD for all pairs
        parfor (cp = 1:nCP, opts.parCores)
            c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
            
            [ld12, ld21] = calc_ld(H, SIG, F, c1, c2);
            
            % fill in output array
            ldArray12a(1,cp,:,w) = ld12; ldArray21a(1,cp,:,w) = ld21;
        end
    end
    
    % Calculate pairwise model
    if ldVersion(2)
        parfor (cp = 1:nCP, opts.parCores)
            c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
            
            % decompose pairwise CPSD
            pairCpsd = thisCpsd([c1,c2],[c1,c2],:);
            [H, SIG] = wilson_sf(pairCpsd, fs);
            
            [ld12, ld21] = calc_ld(H, SIG, F, 1, 2);
            
            % fill in output array
            ldArray12b(1,cp,:,w) = ld12; ldArray21b(1,cp,:,w) = ld21;
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
            ldFeatures12(1, cp,:) = cellfun(@(x) [areaList{c1} '->' areaList{c2} ' ' x], ...
                fStrings, 'UniformOutput', false);
            ldFeatures21(1, cp,:) = cellfun(@(x) [areaList{c2} '->' areaList{c1} ' ' x], ...
                fStrings, 'UniformOutput', false);
        end
    end
    
end
if opts.parCores, delete(pp), end

%% combine values down here
if sum(ldVersion) == 2
    ldArray = {reshape([ldArray12a; ldArray21a], 2*nCP, F, W), ...
        reshape([ldArray12b; ldArray21b], 2*nCP, F, W)};
elseif ldVersion(1)
    ldArray = reshape([ldArray12a; ldArray21a], 2*nCP, F, W);
else
    ldArray = reshape([ldArray12b; ldArray21b], 2*nCP, F, W);
end

ldFeatures = reshape([ldFeatures12; ldFeatures21], 2*nCP, F);
end

function [ldxy, ldyx] = calc_ld(H, SIG, h, x, y)        
Hxy = H(x,y,:);
SIGyx = SIG(y,y)-(SIG(y,x)/SIG(x,x))*SIG(x,y); % partial covariance

Hyx = H(y,x,:);
SIGxy = SIG(x,x)-(SIG(x,y)/SIG(y,y))*SIG(y,x);

ldyx = zeros(1,h); ldxy = zeros(1,h);
for k = 1:h
    ldyx(k) = real(Hxy(:,:,k)*SIGyx*Hxy(:,:, k)');
    ldxy(k) = real(Hyx(:,:,k)*SIGxy*Hyx(:,:, k)');
end
end