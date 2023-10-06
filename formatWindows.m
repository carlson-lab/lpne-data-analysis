function formatWindows(saveFile, useIntervals, centeredWindows, projectFolder, chanFile, fs, ...
windowLength)
% formatWindows
%   Formats data and labels for use in lpne pipeline
%   INPUTS
%   saveFile: name of '.mat' file where you would like to save the
%     formatted data.
%   useIntervals: (optional) boolean indicating whether to only extract
%       data from specific intervals. If 'true', there must be an 'INT_TIME'
%       folder in the project folder containing time interval files for
%       each LFP file.
%   useCenteredWindows: (optional) boolean indicating whether to create
%       centered windows around given timepoints. If 'true', there must be
%       a 'CENTER_TIME' folder in the project folder containing files with
%       a list of millisecond timestamps to create windows around.
%   projectFolder: (optional) where the LFP data is stored. Should have 'Data'
%       and `CHANS` subfolders.
%   chanFile: (optional) filepath of the excel file input containing channel
%       naming information. See `NMFDemo.ipynb`.
%   fs: (optional) sampling rate, in Hz
%   windowLength: (optional) length of largest window size used, in seconds.
%   SAVED VARIABLES
%   data: MxNXP array of the data. M is the #
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

% Fill in default arguments and get input from the user.
if nargin < 2
    useIntervals = false;
end

if nargin < 3
    centeredWindows = false; %default is not centering windows
end

if nargin < 4
    projectFolder = uigetdir('.', 'Select folder containing Data & CHANS subfolders')
    dummy = input(['Make sure areas in channel info file match other ' ...
        'datasets you plan to combine with this one!!\n', 'ENTER to continue']);
end
if nargin < 5
    [chanFile, chanPath] = uigetfile([projectFolder '/*.xls*'], 'Select channel info file');
    chanFile = [chanPath chanFile];
end
if nargin < 6
    inputs = inputdlg({'Enter sampling rate (Hz):'});
    fs = str2double(inputs{1}); % sampling rate, Hz
end
if nargin < 7
    inputs = inputdlg({'Enter the largest window length used in your analysis (s)'});
    windowLength = str2double(inputs{1}); % length of one window (s)
end
pointsPerWindow = fs*windowLength; %the amount of millisecond timepoints in one window

% load channel info (and strip of unwanted ' symbols)
chanData = readtable(chanFile, 'ReadVariableNames', false);
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
if centeredWindows
    centerFolder = [projectFolder '/CENTER_TIME/']; %get folder that contains the information about the centered windows
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

    % check to make sure a valid center file exists. Some contained nans
    % and were omitted
    centerFile = [centerFolder regexprep(filename,'LFP.mat','CENTER.mat')];
    if isfile(centerFile)
        a = 1;
    else
        continue
    end

    if centeredWindows
        % load timsetamp, trial number, and percent progress
        centerFile = [centerFolder regexprep(filename, 'LFP.mat', 'CENTER.mat')];
        load(centerFile)
        timestamps = T;
    end
    

    % extract mouse name and experiment data from filename
    nameParts = split(filename,'_');
    mousename = nameParts{1};
    date = nameParts{2};

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
            elseif centeredWindows
                nWindows = length(timestamps); %the number of timestamps is the number of windows
            else
                nWindows = floor(length(thisChannel)/pointsPerWindow);
            end

            totalWindows = nWindowsParsed + nWindows;
            fileIdx = nWindowsParsed+1:totalWindows;
            thisData = NaN(pointsPerWindow, length(channames), nWindows, 'single');

            labels.allWindows.mouse(fileIdx) = {mousename};
            labels.allWindows.expDate(fileIdx) = {date};
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

            if centeredWindows
                % save centered window labels
                j = 1;
                for i = fileIdx
                    labels.allWindows.timestamps(i) = {timestamps(j)};
                    labels.allWindows.trialNumber(i) = {trialNumber(j)};
                    labels.allWindows.progress(i) = {progress(j)};
                    j = j+1; %we always want j to be from 1:300 while as we go from file to file, i keeps increasing
                end
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
        elseif centeredWindows
            % for each millisecond timepoint, create a centered window
            % around it of the largest window length and concatenate all
            % together
            usableData = zeros(pointsPerWindow, nWindows); %initialize where windows will go
       
            halfWindow = pointsPerWindow/2; %size of half the window
            
            for i = 1:nWindows
                %create all windows for a channel
                thisCentered = thisChannel(timestamps(i)+(1-halfWindow):timestamps(i)+halfWindow); 
                thisCentered = reshape(thisCentered, pointsPerWindow, 1);
                usableData(1:pointsPerWindow, i) = thisCentered;
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




