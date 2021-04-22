function [ldArray, ldFeatures] = linear_directionality_cond(data, areaList, fs, f,...
    opts)
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

opts = fillDefaultOpts(opts);

% generate list of strings for frequencies
fStrings = compose('%d', f)';

% initialize arrays to be filled in loop below
ldArray = zeros(C, C, F, W, 'single');
ldFeatures = cell(C,C, F);
    
% estimate model order (using a subset of the data to conserve memory)
data = permute( data, [2,1,3]);
if W > opts.ordSampleMaxW
    sampledData = data(:,:, randsample(W, opts.ordSampleMaxW));
else
    sampledData = data;
end
[~,BIC] = tsdata_to_infocrit(sampledData, opts.maxOrder, []);
[~,order] = min(BIC);

if opts.parCores, pp = parpool([2 opts.parCores]); end
a = tic;
for w = 1:W
    thisData = data(:, :, w);
    fprintf('Starting window %d: %2.1fs elapsed\n', w, toc(a))
    
    % for each 'sending' region
    for c1 = 1:C
        
        % calculate reduced residuals
        reducedData = thisData([1:(c1-1), (c1+1):C], :);
        [~,~,resid] = tsdata_to_var(reducedData, order);
        
        % generate CPSD for sending region and residuals
        modifiedData = [resid(1:(c1-1),:); thisData(c1, (order+1:end)); resid(c1:end,:)]';
        thisCpsd = cpsd(modifiedData, modifiedData, ...
            opts.window, opts.overlap, f, fs, 'mimo');
        thisCpsd = permute(thisCpsd, [2,3,1]);
        
        % factorize CPSD
        [H, SIG] = wilson_sf(thisCpsd, fs);
        
        % for each 'recieving' region calculate conditional LD
        parfor (c2 = 1:C, opts.parCores)
        %for c2 = 1:C
            if c2==c1, continue, end
            ld12  = calc_ld_cond(H, SIG, F, c2);
            
            % fill in output array
            ldArray(c1,c2,:,w) = ld12;
        end
    end
    
    if w==1
        for c1 = 1:C
            for c2 = 1:C
                % save labels corresponing to each element of feature arrays above
                ldFeatures(c1, c2,:) = cellfun(@(x) [areaList{c1} '->' areaList{c2} ' ' x], ...
                    fStrings, 'UniformOutput', false);
            end
        end
    end
end
if opts.parCores, delete(pp), end
end

function opts = fillDefaultOpts(opts)
if ~isfield(opts,'maxOrder'), opts.maxOrder = 50; end
if ~isfield(opts,'parCores'), opts.parCores = 0; end
if ~isfield(opts,'ordSampleMaxW'), opts.ordSampleMaxW = 1000; end    
end

function ldyx = calc_ld_cond(H, SIG, h, x)
notX = 1:size(SIG,1); notX(x) = [];
Hxy = H(x,notX,:);
SIGyx = SIG(notX,notX)-(SIG(notX,x)/SIG(x,x))*SIG(x,notX); % partial covariance

ldyx = zeros(1,h);
for k = 1:h
    ldyx(k) = real(Hxy(:,:,k)*SIGyx*Hxy(:,:, k)');
end
end
