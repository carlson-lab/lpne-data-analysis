function [pdcArray, dtfArray, pdFeatures] = pdc_dtf(data, areaList, fs, f, opts)
% pdc_dtf
%    Estimates partial directed coherence and directed transfer function
%    from multi-channel timeseries data. This is an estimate of
%    communication from one channel to another.
%
%    INPUTS
%    data: NxAxW array
%        Preprocessed (filtered, averaged, checked for saturation) data.
%         A is the # of areas. N=number of frequency points per
%        window. W=number of time windows.
%    areaList: cell array of strings
%        Ordered list of brain areas from which data was recorded. 
%    fs: scalar
%        Sampling rate (Hz) of the data
%    f: vector
%        Defines the (non-zero) frequencies at which the features should
%        be evaluted.
%    opts: structure
%    FIELDS
%      window: integer or vector
%          If an integer, gives length of Hamming subwindows used in
%          Welch's method for estimating a power spectrum.
%          If a vector, indicates the window that should be used in
%          Welch's method. Can also be left as an empty array '[]' to use
%          the matlab defaule window size.
%      overlap: integer
%          Indicate the overlap between sucessive windows when using
%          Welch's method to estimate a power spectrum. If left as empty
%          (i.e. []) then the matlab default is used.
%      parCores: integer
%          Indicates number of cores to use for parallel computing. If 0,
%          all tasks executed in serial. [Default: 0]

% determine which frequencies to evaluate
f = [0, f];
F = numel(f);

[~,C,W] = size(data);

% generate list of strings for frequencies
fStrings = compose('%d', f)';

% initialize arrays to be filled in loop below
pdcArray = zeros(C, C, F, W, 'single');
dtfArray = zeros(C, C, F, W, 'single');

pdFeatures = cell(C, C, F);

a = tic;
for w = 1:W
    thisData = data(:, :, w);
    fprintf('Starting window %d: %2.1fs elapsed\n', w, toc(a))
    
    thisCpsd = cpsd(thisData, thisData, opts.window, opts.overlap, f, fs, 'mimo');
    thisCpsd = permute(thisCpsd, [2,3,1]);
    
    % Calculate full directionality model
    % decompose full CPSD
    [H, ~] = wilson_sf(thisCpsd, fs);
    
    % calculate A for PDC
    A = zeros(size(H));
    for freq = 1:F
        A(:,:,freq) = inv(H(:,:,freq));
    end
    
    % calculate DTF/PDC for all pairs
    for c = 1:C
        % calculate DTF
        denom_dtf = sum(abs(H(c,:,:)).^2, 2);
        dtfArray(c,:,:, w) = bsxfun(@rdivide, abs(H(c,:,:)), sqrt(denom_dtf));
            
        % calculate PDC
        denom_pdc = sum(abs(A(:,c,:)).^2);
        pdcArray(:,c,:,w) = bsxfun(@rdivide, abs(A(:,c,:)), sqrt(denom_pdc));
    end
    
    if w==1
        for c1 = 1:C
            for c2 = 1:C
                % save labels corresponing to each element of feature arrays above
                pdFeatures(c1, c2,:) = ...
                    cellfun(@(x) [areaList{c1} '->' areaList{c2} ' ' x], ...
                    fStrings, 'UniformOutput', false);
            end
        end
    end
    
end

end