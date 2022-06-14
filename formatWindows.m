function formatWindows(saveFile, useIntervals, startOnly)
% formatWindows
%   Formats data and labels for use in lpne pipeline
%   INPUTS
%   saveFile: name of '.mat' file where you would like to save the
%     formatted data.
%   useIntervals: (optional) boolean indicating whether to only extract
%       data from specific intervals. If 'true' there must be an 'INT_TIME'
%       folder in the project folder containing time interval files for
%       each LFP file.
%   SAVED VARIABLES
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

if nargin < 2
    useIntervals = false;
    if nargin < 3
        startOnly = false;
    end
end

% get inputs
projectFolder = uigetdir('.', 'Select folder containing Data & CHANS subfolders')
dummy = input(['Make sure areas in channel info file match other datasets you plan to ' ...
    'combine with this one!!\n', 'ENTER to continue']);
[chanFile, chanPath] = uigetfile([projectFolder '/*.xls*'], 'Select channel info file');
inputs = inputdlg({'Enter sampling rate (Hz):', 'Enter window length (s)'});
fs = str2double(inputs{1}); % sampling rate, Hz
windowLength = str2double(inputs{2}); % length of one window (s)
pointsPerWindow = fs*windowLength;

% load channel info (and strip of unwanted ' symbols)
chanData = readtable([chanPath chanFile], 'ReadVariableNames', false);
chanData = strrep(chanData{:,:}, "'", '');
channames = chanData(:,1);
chanareas = chanData(:,2);
labels.channel = channames;
labels.channelArea = chanareas;

% for each recording file, load and slice data
chansFolder = [projectFolder '/CHANS/'];
dataFolder = [projectFolder '/Data/'];
dataList = dir([dataFolder '*_LFP.mat']);
nSessions = length(dataList);
if useIntervals
    intFolder = [projectFolder '/INT_TIME/'];
end

% initialize variables to be used in loops below
dataCells = {};
tic
nWindowsParsed = 0;
for k = 1:nSessions
    thisFile = dataList(k);
    if thisFile.isdir, continue, end
    
    % clear channel data from last file and load it from the next file
    clear('-regexp','_\d\d')
    filename = thisFile.name;
    dataPath = [dataFolder filename];
    load(dataPath)
    
    % load channel info; save channel names if needed
    clear('CHANNAMES', 'CHANACTIVE')
    chanFile = [chansFolder regexprep(filename, 'LFP.mat', 'CHANS.mat')];
    load(chanFile, 'CHANNAMES', 'CHANACTIVE')
    
    if useIntervals
        %  load interval info
        intFile = [intFolder regexprep(filename, 'LFP.mat', 'TIME.mat')];
        load(intFile, 'INT_TIME')
        intStart = INT_TIME(1:2:end)*fs + 1;
        intDuration = INT_TIME(2:2:end);
        numIntWindows = floor(intDuration/windowLength);
        intEnd = (intStart -1) + numIntWindows*pointsPerWindow;
    end
    
    % extract mouse name and experiment data from filename
    nameParts = split(filename,'_');
    mousename = nameParts{1};
    date = nameParts{2};
    exp = nameParts{3};
    
    % for every channel, reshape into windows and add to main data
    % array
    channel = who('-regexp','_\d\d');
    C = length(channel);
    for c = 1:C
        channelIdx = strcmp(labels.channel,channel{c});
        thisChannel = eval(channel{c});
        
        % populate labels once for this file
        if c==1
            if useIntervals
                I = length(intStart);
                nWindows = sum(numIntWindows);
            else
                if startOnly
                    nWindows = floor(startOnly*fs/pointsPerWindow);
                else
                    nWindows = floor(length(thisChannel)/pointsPerWindow);
                end
            end
            
            totalWindows = nWindowsParsed + nWindows;
            fileIdx = nWindowsParsed+1:totalWindows;
            thisData = NaN(pointsPerWindow, length(channames), nWindows, 'single');
            
            labels.allWindows.mouse(fileIdx) = {mousename};
            labels.allWindows.expDate(fileIdx) = {date};
            labels.allWindows.exp(fileIdx) = {exp};
            labels.allWindows.time(fileIdx) = 1:nWindows;
            
            if useIntervals
                % save interval labels
                intLabels = zeros(1, nWindows);
                wStart = zeros(1,I); wEnd = zeros(1,I);
                for i = 1:I
                    wStart(i) = sum(numIntWindows(1:(i-1))) + 1;
                    wEnd(i) = sum(numIntWindows(1:i));
                    intLabels(wStart(i):wEnd(i)) = i;
                end
                labels.allWindows.interval(fileIdx) = intLabels;
            end
        end
        
        % skip inactive or unused channels, leaving them as nans
        activeIdx = strcmp(CHANNAMES,channel{c});
        if sum(channelIdx) ~=1 || ~CHANACTIVE(activeIdx)
            continue
        end
        
        if useIntervals
            % extract intervals, slice data for each interval into windows and concatenate
            usableData = zeros(pointsPerWindow, nWindows);
            for i = 1:I
                thisInterval = thisChannel(intStart(i):intEnd(i));
                thisInterval = reshape(thisInterval, pointsPerWindow, numIntWindows(i));
                usableData(:,wStart(i):wEnd(i)) = thisInterval;
            end
        else
            % set whole recording as single interval
            usableData = thisChannel(1:(nWindows*pointsPerWindow));
            usableData = reshape(usableData, pointsPerWindow, nWindows);
        end
        
        thisData(:,channelIdx,:) = usableData;
        
    end
    nWindowsParsed = totalWindows;
    
    dataCells = cat(3,dataCells, {thisData});
    
    fprintf('Day %s of %s loaded. %.1f minutes elapsed\n', ...
        date,mousename,toc/60)
end

% concatenate all files
data = cell2mat(dataCells);

% fill in remaining labels and save
labels.fsRaw = fs;
labels.windowLength = windowLength;

% save data
save(saveFile,'data','labels','-v7.3')
end
