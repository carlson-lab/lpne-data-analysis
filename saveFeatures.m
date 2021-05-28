function saveFeatures(saveFile, options)
% estimate spectral features
%
% INPUTS
% saveFile: name of '.mat' file containing the X and labels
%   variables; and to which the processed data will be saved.
% options: structure of options for running this function
% FIELDS
%   windowOpts: (Optional) a binary vector with length=number of windows,
%     determining which windows to analyze for features.
%   featureList: (Optional) Cell array of strings indicating which
%     features to calculate. The following strings are available options:
%     'power', 'coherence', 'directionality', 'directionality_pairwise',
%     'fft', 'granger'.
%   version: Required if options.featureList is given. Structure contining
%     fields named after each feature listed in feature list, giving the
%     version to be used for that feature (e.g. options.version.power =
%     '1.1'). options.version.directionality corresponds to both
%     directionality and directionality_pairwise features.
%   parCores: (Optional) integer indicating number of cores to use for
%     parallel computing. If 0, all tasks executed in serial.
%   window: integer or vector
%     If an integer, gives length of Hamming subwindows used in
%     Welch's method for estimating a power spectrum.
%     If a vector, indicates the window that should be used in Welch's
%     method. Can also be left as an empty array '[]' to use the
%     matlab defaule window size.
%   overlap: integer
%     Indicate the overlap between sucessive windows when using
%     Welch's method to estimate a power spectrum. If left as empty
%     (i.e. []) then the matlab default is used.
%   mvgcFolder: needed if featureList includes 'granger'. String
%     indicating folder containing MVGC toolbox.
%
% LOADED VARIABLES (from saveFile)
% X: Preprocessed (filtered, averaged, checked for saturation) data. NxAxW array. A is
%       the # of areas. N=number of frequency points per
%       window. W=number of time windows.
% labels: Structure containing labeling infomation for data
%   FIELDS USED
%   fsRaw: sampling frequency of unprocessed data (Hz)
%   area: cell array of labels for each area corresponding to
%       the second dimension of xFft
%   fs: sampling frequency of processed data (Hz)
%   windows: same as allWindows, but with unusable windows eliminated
%   windowLength: length of windows (s)
%
% SAVED VARIABLES
% power: MxNxP matrix of power values where M is frequency, N is brain area,
%     and P is time window
% coherence: MxNxPxQ array to store coherence variables where M is frequency,
%     N is time window, and P and Q are the two brain areas where coherence is
%     calculated
% granger: PxFxW array to store granger causality values. P
%     iterates over directed pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% instant: PxFxW array to store instantaneous causality values. P
%     iterates over undirected pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% directionality: PxFxW array to store 'full' model linear directionality
%     features. P iterates over directed pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% directionality_pairwise: PxFxW array to store pairwise linear directionality
%     features. P iterates over directed pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% fft: fourier transform of X
% labels: See above, with
%   ADDED FIELDS
%   f: integer frequency of processed data
%   powerFeatures: MxN matrix of string labels describing the
%       features represented in power. M is frequency, N is
%       brain area.
%   cohFeatures: MxPxQ array of string labels describing the
%       features represented in coherence. M is frequency, P
%       and Q are the two brain areas where coherence is calculated.
%   gcFeatures: PxF array of string labels describing the
%       features represented in granger. P iterates over
%       directed pairs of regions, F iterates over frequencies.
%   instFeatures: PxF array of string labels describing the
%       features represented in instArray. P iterates over
%       undirected pairs of regions, F iterates over frequencies.
%   ldFeatures: PxF array of string labels describing the
%       features represented in directionality and directionality_pairwise
%       P iterates over directed pairs of regions, F iterates over
%       frequencies.

% original implementation of Welch's method here used 1/4s long windows
ORIG_WELCH_WIN_LEN = 1/4; 

%% load data and prep for feature generation
load(saveFile,'X', 'labels')
fs = labels.fs;

if nargin < 2
    % fill with default parameters
    options=[];
end
options=fillDefaultOpts(options, fs);

% evaluate at every integer frequency up to nyquist
f = 1:floor(fs/2);
nFreq = numel(f);
labels.f = f;

% check if in the format of output from dataSegments algorithm and if so,
% convert to array format as expected
if exist('dataSegments','var')
    X = datautils.segments2array(dataSegments, windowTimes);
elseif exist('X','var')
    % leave data unchanged
else
    error('Data does not appear to be saved as expected in %s\n', saveFile)
end

if ~isempty(options.windowOpts)
    myIdx=find(options.windowOpts==1);
    X=X(:,:,myIdx);
    windowLabels = fieldnames(labels.windows);
    for m = 1:numel(windowLabels)
        labels.windows.(windowLabels{m}) = labels.windows.(windowLabels{m})(myIdx);
    end
end

% convert into 2D matrix NxCW, N is number of windows, CW iterates by brain
% area and then by frequency
[N,C,W] = size(X);
xReshaped = reshape(X, N, C*W);

fStrings = compose('%d', f)';

% check if any features use Welch's method, if so check version. Then
% revert to old Welch's method windowing options if necessary.
welchFeatures = {'power','coherence','directionality'};
calcWFeatures = ismember(welchFeatures, options.featureList);
changeWelchVersion = [false; false; false];
for k = 1:length(calcWFeatures)
   if calcWFeatures(k)
      if ~(str2double(options.version.(welchFeatures{k})(end-2:end)) >= 1.6)
          changeWelchVersion(k) = true;
      end
   end
end
if any(changeWelchVersion)
    warning(['Power, coherence, and directionality versions saveFeatures_1.5 and '...
        'earlier did not use Welch method window options; If you set options.window or'...
        ' options.overlap, they will be overridden.'])
    options.window = round(fs*ORIG_WELCH_WIN_LEN);
    options.overlap = [];
end
    
%% get power spectrum
if any(ismember('power', options.featureList)) && ...
        ~strncmp(options.version.power, 'ld_', 3)

    % Estimate power using Welch's power spectrum estimator
    if strcmp(options.version.power, '1.1')
        scale = 'power';
    else
        scale = 'psd';
    end
    
    power = pwelch(xReshaped, options.window, options.overlap, f,fs, scale);
    % Reshape to MxNxP matrix where M is frequency, N is brain area, and P is time window
    power = reshape(power, [], C, W);
    power = single(abs(power));

    % generate/save matrix where elements name corresponding
    % feature in power matrix
    freqMat = repmat(fStrings, [1 C]);
    areaMat = repmat(labels.area', [nFreq 1]);
    labels.powerFeatures = string(cellfun(@(x,y) [x ' ' y], areaMat, freqMat, ...
        'UniformOutput', false));

    labels.powVersion = options.version.power;

    save(saveFile, 'power', '-append')

elseif any(ismember('power', options.featureList)) && ...
        ~any(ismember('directionality', options.featureList))
    error(['Power cannot be calculated using version %s without adding ''directionality'' to'...
        ' feature list'], options.version.power)
end

%% Get the cross spectra

% Initialize a MxNxPxQ matrix to store coherence variables where M is frequency,
% N is time window, and P and Q are the two brain areas where coherence is
% calculated
if any(ismember('coherence',options.featureList))
    coherence = zeros(nFreq,W,C,C);
    % initialize feature labels to empty strings
    cohFeatures = cell(nFreq,C,C);
    cohFeatures(:) = {''};

    if options.parCores, pp = parpool([2 options.parCores]); end
    a = tic;
    for c1=1:C
        x = squeeze(X(:, c1, :));

        areaList = labels.area;
        parfor (c2=(c1+1):C, options.parCores)
            y = squeeze(X(:, c2, :));

            % use coherence function to calculate coherence between each pair of
            % channels; save to each brain pair location (ie X and Y as well as Y
            % and X)
            cxy = mscohere(x,y, options.window, options.overlap, f,fs);
            coherence(:,:,c2,c1) = cxy; % symmetric

            % save feature labels for this pair of regions
            cohFeatures(:,c2,c1) = cellfun(@(x) [areaList{c1} '-' areaList{c2} ' ' x], fStrings, ...
                'UniformOutput', false);
        end

        fprintf('Evaluated remaining %s coherence(s): %2.1fs elapsed\n', areaList{c1}, toc(a))
    end
    if options.parCores, delete(pp), end

    labels.cohFeatures = string(cohFeatures);

    labels.cohVersion = options.version.coherence;

    coherence=single(abs(coherence));
    % save coherence results to data
    save(saveFile, 'coherence', '-append')
end

%% Get Granger causality features

if any(ismember('granger', options.featureList))
    mvgcStartupScript = [options.mvgcFolder '/startup.m'];
    run(mvgcStartupScript)

    % only high pass filter data for versions <1.5
    if strcmp(options.version.granger, 'saveFeatures_1.5')
        X_filt = double(X);
    else
        d = designfilt('highpassiir', 'PassbandFrequency', 1/ORIG_WELCH_WIN_LEN, ...
            'StopbandFrequency', 1/labels.windowLength , 'SampleRate', fs);
        X_filt = filtfilt(d, double(X));
    end

    % generate Granger causality values matrix in the form MxPxQ, where M
    % iterates over pairs of brain regions, P is frequency, and Q is time
    % window.
    [granger, gcFeatures, instant, instFeatures] = g_causality(X_filt, labels.area, fs, ...
                                                               options);
    granger = single(granger);
    labels.gcFeatures = string(gcFeatures);
    instant = single(instant);
    labels.instFeatures = string(instFeatures);

    labels.gcVersion = options.version.granger;

    save(saveFile, 'granger', 'instant', '-append')
end

%% Get linear directionality features
% Check if full or pairwise directionality features are to be calculated.
directionFeatures = ismember({'directionality', 'directionality_pairwise'}, ...
    options.featureList);
if any(directionFeatures)
    X = double(X);
    
    % generate linear directionality values matrix in the form MxPxQ, where M
    % iterates over pairs of brain regions, P is frequency, and Q is time
    % window.
    [directionality, ldFeatures, S] = linear_directionality(X, labels.area, fs, f,...
        directionFeatures, options);
    
    labels.ldFeatures = string(ldFeatures);
    labels.ldVersion = options.version.directionality;
    
    % Save calculated features
    if sum(directionFeatures) == 2
        directionality_pairwise = single(directionality{2});
        directionality = single(directionality{1});
        save(saveFile, 'directionality', 'directionality_pairwise', '-append')
    elseif directionFeatures(1)
        directionality = single(directionality);
        save(saveFile, 'directionality','-append')
    else
        directionality_pairwise = single(directionality);
        save(saveFile, 'directionality_pairwise','-append')
    end  
    
    % Check if power values should be saved from the CPSD generated in ld
    % calculations
    if any(ismember('power', options.featureList)) && ...
            strncmp(options.version.power, 'ld_saveFeatures', 15)
        power = zeros(nFreq, C, W);
        for k =1:C
            power(:,k,:) = S(k,k,:,:);
        end

        % generate/save matrix where elements name corresponding
        % feature in power matrix
        freqMat = repmat(fStrings, [1 C]);
        areaMat = repmat(labels.area', [nFreq 1]);
        labels.powerFeatures = string(cellfun(@(x,y) [x ' ' y], areaMat, freqMat, ...
        'UniformOutput', false));

        labels.powVersion = options.version.power;

        save(saveFile, 'power','-append')
    end
end

if ismember('psi', options.featureList)
    segleng = fs; % 1 Hz frequency resolution
    fBins = [f(1:end-1)', (f(1:end-1)'+1), (f(1:end-1)'+2)];
    
    % create feature labels
    psiFeatures = cell(C,C,nFreq-1);
    psiFeatures(:) = {''};
    areaList = labels.area;
    for c1 = 1:C
        for c2 = 1:C
            if c1==c2, continue, end
            % save feature labels for this pair of regions
            psiFeatures(c1,c2,:) = cellfun(@(x) [areaList{c1} '~>' areaList{c2} ' ' x], ...
                fStrings(1:end-1), 'UniformOutput', false);
        end
    end
    
    labels.psiFeatures = string(psiFeatures);
    labels.psiVersion = 'saveFeatures_1.6';
    
    if options.parCores, pp = parpool([2 options.parCores]); end
    
    % calculate psi for eac window
    psi = zeros(C,C,nFreq-1,W, 'single');
    parfor (w = 1:W, options.parCores)
        thisData = X(:,:,w);
        psi(:,:,:,w) = data2psi(thisData, segleng, [], fBins);
    end
    if options.parCores, delete(pp), end
    
    save(saveFile, 'psi', '-append')
end

if any(ismember({'pdc','dtf'}, options.featureList))
    % generate linear directionality values matrix in the form MxPxQ, where M
    % iterates over pairs of brain regions, P is frequency, and Q is time
    % window.
    X = double(X);
    
    [pdc, dtf, pdFeatures] = pdc_dtf(X, labels.area, fs, f, options);
    
    labels.pdFeatures = string(pdFeatures);
    labels.pdVersion = 'saveFeatures_1.6';
    
    save(saveFile, 'pdc', 'dtf', '-append')
end


%% Take Fourier transform of data
if any(ismember('fft', options.featureList))
    Ns = ceil(N/2);
    if strcmp(options.version.fft, '1.1')
        scale = 1/sqrt(N);
    else
        scale = 1/N;
    end
    xFft = scale*fft(double(X));
    xFft = 2*(xFft(2:Ns+1,:,:));

    labels.fftVersion = options.version.fft;

    save(saveFile,'xFft','-append');
end

%% Save labels to JSON
save(saveFile,'labels','-append');
datautils.saveJson(saveFile, labels)

end

function opts = fillDefaultOpts(opts, fs)
    if ~isfield(opts,'windowOpts'), opts.windowOpts = []; end
    if ~isfield(opts,'featureList')
        opts.featureList = {'power','coherence','directionality'};
        opts.version.power = 'saveFeatures_1.6';
        opts.version.coherence = 'saveFeatures_1.6';
        opts.version.directionality = 'saveFeatures_1.6';
    end
    if ~isfield(opts,'parCores'), opts.parCores = 0; end
    if ~isfield(opts, 'window'), opts.window = rectwin(round(fs*2/5)); end
    if ~isfield(opts, 'overlap'), opts.overlap = []; end
    if ~isfield(opts,'mvgcFolder'), opts.mvgcFolder = '~/lpne-data-analysis/mvgc'; end
end
