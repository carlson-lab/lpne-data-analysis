function saveFeatures(saveFile, options)
% estimate spectral features
% INPUTS
% saveFile: name of '.mat' file containing the data and labels
%   variables. And to which the processed data will be saved
% options: structure of options for running this function
% FIELDS
%   windowOpts: (Optional) a binary vector with length=number of windows,
%     determining which windows to analyze for features.
%   featureList: (Optional) Cell array of strings indicating which
%     features to calculate
%   parCores: (Optional) integer indicating number of cores to use for
%     parallel computing. If 0, all tasks executed in serial.
%   mvgcFolder: needed if featureList includes 'granger'. String
%     indicating folder containing MVGC toolbox.
%
% LOADED VARIABLES
% labels: Structure containing labeling infomation for data
%   FIELDS USED
%   fsRaw: sampling frequency of unprocessed data (Hz)
%   area: cell array of labels for each area corresponding to
%       the second dimension of xFft
%   fs: sampling frequency of processed data (Hz)
%   windows: same as allWindows, but with unusable windows eliminated
%   X: Preprocessed (filtered, averaged, checked for saturation) data. NxAxW array. A is
%       the # of areas. N=number of frequency points per
%       window. W=number of time windows.
%   windowLength: length of windows (s)
% dataSegments (optional): output data from preprocessing through
%   dataSegments method.
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
% labels: See above, with
%   ADDED FIELDS
%   f: integer frequency of processed data
%   powerFeatures: MxN matrix of string labels describing the
%       features represented in labels.power. M is frequency, N is
%       brain area.
%   cohFeatures: MxPxQ array of string labels describing the
%       features represented in labels.coherence. M is frequency, P
%       and Q are the two brain areas where coherence is calculated.
%   gcFeatures: PxF array of string labels describing the
%       features represented in labels.granger. P iterates over
%       directed pairs of regions, F iterates over frequencies.
%   instFeatures: PxF array of string labels describing the
%       features represented in instArray. P iterates over
%       undirected pairs of regions, F iterates over frequencies.
FEATURES_VERSION = '1.1';
WELCH_WIN_LEN = 1/4; % quarter of a second (frequency resolution of 4Hz)
    
if nargin < 2
    % fill with default parameters
    options=[];
end
options=fillDefaultOpts(options);

fprintf('Ignore the following warning(s) \n')
load(saveFile,'X','dataSegments','labels', 'windowTimes')
fs = labels.fs;

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

[N,C,W] = size(X);
xReshaped = reshape(X, N, C*W);
windowSize = round(fs*WELCH_WIN_LEN); 

fStrings = compose('%d', f)';
%% get power spectrum
if any(ismember('power',options.featureList))
    % convert into 2D matrix NxCW, N is number of windows, CW iterates by brain
    % area and then by frequency

    % Estimate power using Welch's power spectrum estimator
    power = pwelch(xReshaped, windowSize,[], f,fs, 'power');
    % Reshape to MxNxP matrix where M is frequency, N is brain area, and P is time window 
    power = reshape(power, [], C, W);
    power = single(abs(power));
    
    % generate/save matrix where elements name corresponding
    % feature in power matrix  
    freqMat = repmat(fStrings, [1 C]);
    areaMat = repmat(labels.area', [nFreq 1]);
    labels.powerFeatures = string(cellfun(@(x,y) [x ' ' y], areaMat, freqMat, ...
        'UniformOutput', false));
    
    labels.powVersion = FEATURES_VERSION;
    
    save(saveFile, 'power', '-append')
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
            cxy = mscohere(x,y, windowSize,[], f,fs);
            coherence(:,:,c2,c1) = cxy; % symmetric
            
            % save feature labels for this pair of regions
            cohFeatures(:,c2,c1) = cellfun(@(x) [areaList{c1} '-' areaList{c2} ' ' x], fStrings, ...
                'UniformOutput', false);
        end

        fprintf('Evaluated remaining %s coherence(s): %2.1fs elapsed\n', areaList{c1}, toc(a))
    end
    if options.parCores, delete(pp), end
    
    labels.cohFeatures = string(cohFeatures);
    
    labels.cohVersion = FEATURES_VERSION;
    
    coherence=single(abs(coherence));
    % save coherence results to data
    save(saveFile, 'coherence', '-append')
end

%% Get Granger causality features

if any(ismember('granger', options.featureList))
    mvgcStartupScript = [options.mvgcFolder '/startup.m'];
    run(mvgcStartupScript)
    
    % generate Granger causality values matrix in the form MxPxQ, where M 
    % iterates over pairs of brain regions, P is frequency, and Q is time
    % window.
    d = designfilt('highpassiir', 'PassbandFrequency', 1/WELCH_WIN_LEN, ...
                   'StopbandFrequency', 1/labels.windowLength , 'SampleRate', fs);
    X_filt = filtfilt(d, double(X));
    [granger, gcFeatures, instant, instFeatures] = g_causality(X_filt, labels.area, fs, ...
                                                               options);
    granger = single(granger);
    labels.gcFeatures = string(gcFeatures);
    instant = single(instant);
    labels.instFeatures = string(instFeatures);
    
    labels.gcVersion = FEATURES_VERSION;
    
    save(saveFile, 'granger', 'instant', '-append') 
end

%% Take Fourier transform of data
if any(ismember('fourier', options.featureList))
    Ns = ceil(N/2);
    xFft = 1/sqrt(N)*fft(double(X));
    xFft = 2*(xFft(2:Ns+1,:,:));
    
    labels.fftVersion = FEATURES_VERSION;
    
    save(saveFile,'xFft','-append');
end

%% Save labels to JSON
save(saveFile,'labels','-append');
datautils.saveJson(saveFile, labels)

end

function opts = fillDefaultOpts(opts)
    if ~isfield(opts,'windowOpts'), opts.windowOpts = []; end
    if ~isfield(opts,'featureList')
        opts.featureList = {'power','coherence','granger'};
    end
    if ~isfield(opts,'parCores'), opts.parCores = 0; end
    if ~isfield(opts,'mvgcFolder'), opts.mvgcFolder = '~/lpne-data-analysis/mvgc'; end
end
