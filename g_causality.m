function [gcArray, gcFeatures] = g_causality(data, areaList, fs, options)
% Frequency-domain granger causality
%    Estimates frequency-domain granger causality from x->y and y->x.
%    Estimates are given for all frequency bands given by f.
%
%    INPUTS
%    data: Preprocessed (filtered, averaged, checked for saturation) data.
%        NxAxW array. A is the # of areas. N=number of frequency points
%        per window. W=number of time windows.
%    fs: Sampling rate (Hz)
%    areaList: ordered list of brain areas from which data was recorded
%    options: structure of options for this function. If any not included,
%        will be filled in with default values.
%    FIELDS
%      separateWindows: boolean indicating whether to model
%          windows seperately [Default: true]
%      maxOrder: maximum model order for model order estimation [Defaut: 50]
%      parCores: integer indicating number of cores to use for%
%          parallel computing. If 0, all tasks executed in serial. [Default: 0]
%      ordSampleMaxW: maximum number of windows to use in estimating the
%          order of the AR models used to calculate Granger values.
%    OUTPUTS
%    gcArray: spectral granger causality estimates
%    gcFeatures: PxF array of string labels describing the
%        features represented in gcArray. P iterates over
%        directed pairs of regions, F iterates over frequencies.
%    instArray: instantaneous causality estimates
%    instFeatures: PxF array of string labels describing the
%        features represented in instArray. P iterates over
%        undirected pairs of regions, F iterates over frequencies.
%
% This function requires the MVGC toolbox. Portions of this code were taken
% from mvgc_demo script of that toolbox.
% Toolbox can be accessed at: https://users.sussex.ac.uk/~lionelb/MVGC/

if nargin < 2
    % fill with default parameters
    options=[];
end
options=fillDefaultOpts(options);

data = permute( data, [2,1,3]);

% determine which frequencies to evaluate
F = floor(fs/2) + 1;
f = sfreqs(F-1, fs);

if options.separateWindows
    % loop over each window seperately
    [C,~,W] = size(data);
else
    % run through loop once using all windows
    W = 1;
    C = size(data,1);
end

% number of directed pairs
nCP = 2*nchoosek(C, 2)';

% initialize arrays to be filled in loop below
gcArray = zeros(nCP, F, W, 'single');

% estimate model order (using a subset of the data to conserve memory)
if W > options.ordSampleMaxW
    sampledData = data(:,:, randsample(W, options.ordSampleMaxW));
else
    sampledData = data;
end
[~,~,~,order] = tsdata_to_infocrit(sampledData, options.maxOrder, [], false);

if options.parCores, pp = parpool([2 options.parCores]); end
a = tic;
parfor (w = 1:W, options.parCores)
%for w = 1:W
    if options.separateWindows
        thisData = data(:, :, w);
        %fprintf('Starting window %d: %2.1fs elapsed\n', w, toc(a))
    else
        error('Not built for all windows together!')
    end

    % estimate VAR model
    [A,SIG] = tsdata_to_var(thisData, order);
    if warn_if(isbad(A),'in full regression - regression failed'), end

    % estimate autocovariance from VAR model
    [G, info] = var_to_autocov(A, SIG, options.acMaxLags, options.acDecayTol);
    %var_info(info, true);
    if info.error
        fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
        continue
    end
    if info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
        fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)',info.acminlags,realpow(info.rho,info.aclags));
    end

    thisGC = autocov_to_spwcgc(G, F-1);
    assert(isreal(thisGC))
    gcArray(:,:,w) = reshape(thisGC(~isnan(thisGC)), [], F);
    
    fprintf('window completed - %.3fs elapsed\n', toc(a))
end
if options.parCores, delete(pp), end

toc(a)

% save labels corresponing to each element of feature arrays above
% generate list of strings for frequencies
fStrings = reshape(string(compose('%d', f)), 1, 1, []);
preFeatures = areaList + "->" + areaList' + " " + fStrings;
diagMask = repmat(eye(length(areaList), 'logical'), 1, 1, length(fStrings));
preFeatures(diagMask) = "";
gcFeatures = reshape(preFeatures(preFeatures~=""), [], F);
end

function opts = fillDefaultOpts(opts)
if ~isfield(opts,'maxOrder'), opts.maxOrder = 20; end
if ~isfield(opts,'separateWindows'), opts.separateWindows = true; end
if ~isfield(opts,'parCores'), opts.parCores = 0; end
if ~isfield(opts,'ordSampleMaxW'), opts.ordSampleMaxW = 1000; end
if ~isfield(opts,'acMaxLags'), opts.acMaxLags = 3000; end
if ~isfield(opts,'acDecayTol'), opts.acDecayTol = 1e-5; end
end
