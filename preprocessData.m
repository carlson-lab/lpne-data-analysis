function preprocessData(saveFile,dataOpts)
% preprocessData
%   Preprocesses data saved in saveFile.
%   INPUTS
%   saveFile: name of '.mat' file containing the data and labels
%       variables. And to which the fully preprocessed data will be saved
%   dataOpts: optional input. see description in saved variables section. You
%         can  specify some or none of the fields for this structure, and the
%      rest will be filled with default values.
%   LOADED VARIABLES
%   data: MxNXP array of the data for each delay length. M is the #
%       of time points. N is the # of channels. P is the # of
%       windows. All elements corresponding to data that was not
%       saved (i.e. missing channel) should be marked with NaNs.
%   labels: Structure containing labeling infomation for data
%       FIELDS
%       channel: cell arrays of names for channels
%       channelArea: cell array w/ same size as channel giving area
%           assignment for each channel
%       fsRaw: sampling frequency of unprocessed data (Hz)
%       windowLength: length of data windows in seconds
%       allWindows: structure containing relevant labels pertaining to
%           individual windows. Each field should be a vector / array with
%           one element corresponding to each window. Suggested fields:
%           date, etc. Must contain 'mouse', 'expDate', and 'time' fields.
%   SAVED VARIABLES
%   labels: see above
%       ADDED FIELDS
%       area: cell array of labels for each area corresponding to
%           the second dimension of xFft
%       fs: sampling frequency of processed data (Hz)
%       s: frequency space labels of fourier transformed data
%       windows: same as allWindows, but with unusable windows eliminated
%   saturatedPoint: cell vector, where each cell contains logical MxCxW array
%       of the data for each delay length indicating whether a data point is
%       considered saturated (i.e. there is a signal artifact). M=points
%       per time window. C=number of channels. W=number of time windows.
%   averagedData: MxAxW array of the data. Contains average signal
%       for each area. A is the # of areas. M=points
%       per time window. W=number of time windows.
%   (optionally) Xfft: Fourier transform of the data
%   dataOpts: Data preprocessing options.
%       FIELDS
%       subSampFact: subsampling factor
%       normWindows: String indication pf how to normalize
%           individual windows. If 'window', each window is normalized. If 'whole',
%           dataset is normalized as a whole, and individual windows are
%           still mean subtracted. Can also be description 'mouse' to
%           normalize by channel for each mouse, or 'day' to
%           normalize by day for each channel for each mouse. Can
%           be set to 'none' if no normalization is desired. If no
%           option is recognized, the data will be normalized by window.
%           Default is 'day'.
%       transform: whether or not to take and save the Fourier transform of the data
%       satThresh: method for choosing saturation detection threshold.can be a
%           numerical vector or a string. If a vector, each
%           element gives the saturation detection threshold for the
%           corresponding area in the data. If a string, uses a method for
%           automatic detection of thresholds: 'MAD' indicates to used
%           median absolute deviation method to detect outlier points;
%           'SD' indicates to use standard deviation.
PREPROCESS_VERSION = '1.0';

    % exclude saturated data and condense multiple channels from the same brain
    % region to one signal
    load(saveFile, 'data', 'labels');
    data = single(data); % convert to single for memory efficiency
    
    % make sure necessary information has been saved
    assert(isfield(labels.allWindows, 'mouse'), ['Please resave data with mouse identity ' ...
                        'labels in \labels.allWindows.mouse'])
    assert(isfield(labels.allWindows, 'time'), ['Please resave data with time ' ...
                        'labels in \labels.allWindows.time'])

    if nargin < 2
        % fill with default parameters
        dataOpts = [];
    end
    dataOpts = fillDefaultDopts(dataOpts, labels);
    
    [saturatedPoint, dataOpts] = computeSaturation(data, labels, dataOpts);
    [X,labels] = averageAreas(labels,data,saturatedPoint);
    clear('data')
    
    %initialize variables
    fs = labels.fsRaw;
    sampsPerWindow = labels.windowLength * fs;
    nAreas = size(X,2);
    
    % subsample data using subsampling frequency
    % make sure subsample factor is a factor of sampling frequency
    if mod(fs,dataOpts.subSampFact)
        error('Subsampling factor is not a factor of sampling frequency.')
    end

    ptsPerWindow = floor(sampsPerWindow/dataOpts.subSampFact);

    fsFinal = fs/dataOpts.subSampFact;

    %% Preprocessing

    % 60Hz notch filters (+harmonics)
    % stop frequencies at multiples of 60 Hz +/- 2.5 Hz, allow to pass
    % anything greater than 1 Hz away.

    Rp=0.5; % set maximum passband distortion
    Rs=30; % set minimum stopband attenuation
    X = double(X);
    for f = 60:60:floor(fsFinal/2)
        Wp = [(f-5) (f+5)]*2/fs; % passband
        Ws = [(f-2.5) (f+2.5)]*2/fs; % stopband
        [n,Wn] = buttord(Wp,Ws,Rp,Rs);
        [z,p,k] =  butter(n,Wn,'stop');
        [sos,g] = zp2sos(z,p,k);
        X = filtfilt(sos,g,X);
    end

    % subsample (using custom function based off 'decimate')
    X = datautils.decimateiir(X,dataOpts.subSampFact);

    % normalize data based on user-defined option .normWindows, options include
    % 'window','whole', 'mouse', 'day'

    switch dataOpts.normWindows
      case 'window'
        % normalize by windows
        X = zscore(X);

      case 'none'
        X = X;

      case 'whole'
        % zero mean each window
        X = bsxfun(@minus, X, mean(X));      

        % divide by sd to get unit variance
        X_all_w = reshape(permute(X, [1 3 2]), [], size(X,2));
        X_std = std(X_all_w);
        X = bsxfun(@rdivide, X, X_std); 
        dataOpts.normFact = X_std;
        
      case 'mouse'
        % normalize by channel for each mouse
        mouse = unique(labels.allWindows.mouse);
        for m = mouse(:)'
            thisIdx = ismember(labels.windows.mouse,m);
            xThisMouse = X(:,:,thisIdx);
            % windows still have zero mean
            xZeroMean = bsxfun(@minus,xThisMouse,mean(xThisMouse));
            % set Std Dev of each channel to 1
            xZeroMean = reshape(permute(xZeroMean,[1,3,2]),[],nAreas);
            xNorm = bsxfun(@rdivide,xZeroMean,std(xZeroMean));
            X(:,:,thisIdx) = permute(reshape(xNorm,ptsPerWindow,[], ...
                                             nAreas),[1,3,2]);
        end
      case 'day'
        % normalize by channel for each mouse
        mouse = unique(labels.allWindows.mouse);
        day = unique(labels.allWindows.expDate);
        for m = mouse(:)'
            for d = day(:)'
                thisIdx = ismember(labels.windows.mouse,m) & ...
                          ismember(labels.windows.expDate,d);
                if ~sum(thisIdx), continue, end
                xThisDay = X(:,:,thisIdx);
                % windows still have zero mean
                xZeroMean = bsxfun(@minus,xThisDay,mean(xThisDay));
                % set Std Dev of each channel to 1
                xZeroMean = reshape(permute(xZeroMean,[1,3,2]),[],nAreas);
                xNorm = bsxfun(@rdivide,xZeroMean,std(xZeroMean));
                X(:,:,thisIdx) = permute(reshape(xNorm,ptsPerWindow,[], ...
                                                 nAreas),[1,3,2]);
            end
        end
      otherwise
        warning('data.normWindows value not recognized. Normalizing by window')
        X = zscore(X);
    end


    % Save frequency labels
    labels.s = (fsFinal/ptsPerWindow):(fsFinal/ptsPerWindow):ceil(fsFinal/2);
    labels.fs = fsFinal;

    % if transform is wanted, take Fourier transform and save to file
    if dataOpts.transform
        % take Fourier transform
        Ns = ceil(ptsPerWindow/2);
        xFft = 1/sqrt(ptsPerWindow)*fft(X);
        xFft = 2*(xFft(2:Ns+1,:,:));
        save(saveFile,'xFft','-append');
    end

    % save final data to file
    X = single(X);
    labels.preprocessVersion = PREPROCESS_VERSION;
    save(saveFile,'dataOpts','labels','X','-append')
end


function dataOpts=fillDefaultDopts(dataOpts, labels)
%fill in dataOpts with default data options:
    if ~isfield(dataOpts,'subSampFact'), dataOpts.subSampFact = 2; end 
    if ~isfield(dataOpts,'normWindows'), dataOpts.normWindows = 'whole'; end
    if ~isfield(dataOpts,'transform'), dataOpts.transform = 0; end
    if ~isfield(dataOpts,'satThresh'), dataOpts.satThresh = 'MAD'; end
    
    dataOpts.windowLength = labels.windowLength;
    dataOpts.fsRaw = labels.fsRaw;
end

function [averagedData, labels] = averageAreas(labels,data,saturatedPoint)
%AverageChannels - Compute average signal from channels for each area
%   Computes the average signal for each area at
%   all time points in the data matrix by combining channels in that area.
%   Any points marked as saturated in each channel are excluded.
%   INPUTS
%   data: cell vector, where each cell contains MxNXP array of the data for each
%       delay length. M is the # of time points. N is the # of channels. P
%       is the # of windows
%   labels: Structure containing labeling infomation for data
%   saturatedPoints: cell vector, where each cell contains locical MxNXP array
%       of the data for each delay. length indicating whether a data point is
%       considered saturated (i.e. there is a signal artifact).  M is the # of
%       time points. N is the # of channels. P is the # of windows
%   OUTPUTS
%   averagedData: cell vector, where each cell contains MxAxP array of the data
%       for each delay. Each array contains average signal for each area. A is
%       the # of areas. If all channels for an area are saturated at a certain
%       point, value at that element is set to NaN.
%   labels: Structure containing labeling infomation for data

% Get names for each of the unique recording areas
    area = unique(labels.channelArea);
    %% If there are any channels that need to be removed for any reason, do
    % that here in remove and then uncomment the following lines:

    % remove = ["Cg_Cx_L","IL_Cx_L","PrL_Cx_L","IL_Cx_R","PrL_Cx_R"];
    % conts = cellfun(@(c)ismember(c,remove),area,'UniformOutput',false);
    % idx = []; 
    % for i=1:length(area) 
    %     if conts{i}==1 
    %         idx=[idx i];
    %     end, 
    % end
    % for i=length(idx):-1:1
    % area(idx(i))=[];
    % end

    %% get sizes for averagedData array
    nAreas = numel(area);
    nPts = size(data,1);
    nWin = size(data,3);

    % initialize averagedData array with zeros
    averagedData = zeros([nPts,nAreas,nWin],'single');

    % replace each saturated point with nan
    data(saturatedPoint) = nan;
    bad={};
    % compute averaged data for each area, ignoring nan
    for a = 1:nAreas
        chanInArea = strcmp(labels.channelArea,area{a});
        areaData = data(:,chanInArea,:);
        averagedData(:,a,:) = mean(areaData,2,'omitnan');
        
        if sum(sum(isnan(averagedData(:,a,:))))/numel(averagedData(:,a,:))>0.75
            bad{end+1}=area{a};
        end
    end
    if ~isempty(bad)
        for i=1:numel(bad)
            error("Won't work, too many unusable points in %s, stopping", bad{i})
        end
        
    end
    % remove windows w/ NaNs
    nansInData = isnan(averagedData);
    winHasNan = logical(squeeze(sum(sum(nansInData,2),1)));
    averagedData(:,:,winHasNan) = [];
    windowLabels = fieldnames(labels.allWindows);
    labels.windows = labels.allWindows;
    for f = 1:numel(windowLabels)
        labels.windows.(windowLabels{f})(winHasNan) = [];
    end

    labels.area = area;
end


