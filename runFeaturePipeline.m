function runFeaturePipeline(saveFile, dataOpts, featureOpts, trainFile)
% run all portions of pipeline required to generate frequency based features from individual recording files
%   INPUTS
%   saveFile: name of '.mat' file where you would like to save the
%     formatted data.
%   dataOpts: optional input. You
%      can  specify some or none of the fields for this structure, and the
%      rest will be filled with default values.
%      FIELDS
%       subSampFact: subsampling factor
%       normWindows: String indication pf how to normalize
%           individual windows. If 'window', each window is normalized. If 'whole',
%           dataset is normalized as a whole, and individual windows are
%           still mean subtracted. Can also be description 'mouse' to
%           normalize by channel for each mouse, or 'day' to
%           normalize by day for each channel for each mouse. Can
%           be set to 'none' if no normalization is desired. If no
%           option is recognized, the data will be normalized by window.
%           Default is 'whole'.
%       transform: whether or not to take and save the Fourier transform of the data
%       satThresh: method for choosing saturation detection threshold.can be a
%           numerical vector or a string. If a vector, each
%           element gives the saturation detection threshold for the
%           corresponding area in the data. If a string, uses a method for
%           automatic detection of thresholds: 'MAD' indicates to used
%           median absolute deviation method to detect outlier points;
%           'SD' indicates to use standard deviation.
%   featureOpts: structure of options for running this function. You
%      can  specify some or none of the fields for this structure, and the
%      rest will be filled with default values.
%     FIELDS
%       windowOpts: (Optional) a binary vector with length=number of windows,
%         determining which windows to analyze for features.
%       featureList: (Optional) Cell array of strings indicating which
%         features to calculate
%       parCores: (Optional) integer indicating number of cores to use for
%         parallel computing. If 0, all tasks executed in serial.
%       mvgcFolder: needed if featureList includes 'granger'. String
%         indicating folder containing MVGC toolbox.
%   trainFile (optional): indicates location of saved/preprocessed training data if
%       this data will be backprojected into a model trained on other data.
%       Do not set if this data is for training. If set, uses the training
%       file to determine the saturation thresholds and normalization
%       constants for preprocessing.
  
if nargin < 4
    trainFile = false;
    input(['No training file indicated. If you plan to use this data for backprojection, '...
        'you might want to set the trainFile parameter.\n'...
        'If the indicated data is for training a new model, make sure you have '...
        'chosen the correct normalization method.\nWhich type of normalization to use depends on '...
        'the experimental question; If you don''t know why this is, ask someone who does.\n'...
        'See runFeaturePipeline documentation and preprocessData documentation for more details.\n'...
        'Press ENTER to continue'])
    if nargin < 3
        featureOpts = []
        if nargin < 2
            dataOpts = []
        end
    end
end

formatWindows(saveFile)

if trainFile
    train = load(trainFile, 'dataOpts');
    dataOpts.satThresh = train.dataOpts.satThresh;
    dataOpts.normWindows = train.dataOpts.normWindows;
    dataOpts.subSampFact = train.dataOpts.subSampFact;
    if isfield(train.dataOpts, 'normFact')
        dataOpts.normFact = train.dataOpts.normFact;
    end
end
  
preprocessData(saveFile, dataOpts)

saveFeatures(saveFile, featureOpts)
