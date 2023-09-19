function formatWindows(saveFile, useIntervals, projectFolder, chanFile, fs, windowLength)

% Fill in default arguments and get input from the user.
if nargin < 2
    useIntervals = false;
end
if nargin < 3
    projectFolder = uigetdir('.', 'Select folder containing Data & CHANS subfolders');
    dummy = input(['Make sure areas in channel info file match other ' ...
        'datasets you plan to combine with this one!!\n', 'ENTER to continue']);
end
if nargin < 4
    [chanFile, chanPath] = uigetfile([projectFolder '/*.xls*'], 'Select channel info file');
    chanFile = [chanPath chanFile];
end
if nargin < 5
    inputs = inputdlg({'Enter sampling rate (Hz):'});
    fs = str2double(inputs{1}); % sampling rate, Hz
end
if nargin < 6
    inputs = inputdlg({'Enter window length (s)'});
    windowLength = str2double(inputs{1}); % length of one window (s)
end
pointsPerWindow = fs*windowLength;

% Load channel info (and strip of unwanted ' symbols)
chanData = readtable(chanFile, 'ReadVariableNames', false);
chanData = strrep(chanData{:,:}, "'", '');
channames = chanData(:,1);
chanareas = chanData(:,2);
labels.channel = channames;
labels.channelArea = chanareas;

% For each recording file, load and slice data
chansFolder = [projectFolder '/CHANS/'];
dataFolder = [projectFolder '/Data/'];
dataList = dir([dataFolder '*_LFP.mat']);
nSessions = length(dataList);
if useIntervals
    intFolder = [projectFolder '/INT_TIME/'];
end

% Initialize variables to be used in loops below
dataCells = {};
tic
nWindowsParsed = 0;

totalWindows = 0; %this line was added and is not in the original

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

    % Load the frame file for the current LFP file, ALL OF THIS IS NEW
    frameFile = [projectFolder '/Frame/' regexprep(thisFile.name, 'LFP.mat', 'FRAME.mat')];
    if exist(frameFile, 'file')
        load(frameFile, 'frames')
    else
        error('Frame file not found for %s', thisFile.name);
    end

    % Calculate the window boundaries for each frame, STILL NEW
    intStartBefore = frames - windowLength * fs;
    intEndBefore = frames - 1;
    intStartAfter = frames + 1;
    intEndAfter = frames + windowLength * fs;
    
    % Remove any negative intStart values and corresponding intEnd values,
    % STILL NEW
    invalidIndices = intStartBefore <= 0;
    intStartBefore(invalidIndices) = [];
    intEndBefore(invalidIndices) = [];
    intStartAfter(invalidIndices) = [];
    intEndAfter(invalidIndices) = [];

    nWindows = 2 * length(intStartBefore); % NEW

    % Extract mouse name and experiment date from filename
    nameParts = split(thisFile.name,'_');
    mousename = nameParts{1};
    date = nameParts{2};

    % Create data structure for this file's windows, THIS IS NEW
    thisData = NaN(pointsPerWindow, length(channames), nWindows);

    % For every channel, reshape into windows and add to main data array
    channel = who('-regexp','_\d\d');
    C = length(channel);
    for c = 1:C
        channelIdx = strcmp(labels.channel, channel{c});
        thisChannel = eval(channel{c});

        % Skip inactive or unused channels, leaving them as NaNs
        activeIdx = strcmp(channames, channel{c});
        if sum(channelIdx) ~= 1
            continue;
        end

        usableDataBefore = zeros(pointsPerWindow, nWindows/2);
        usableDataAfter = zeros(pointsPerWindow, nWindows/2);

        for i = 1:nWindows/2
            thisIntervalBefore = thisChannel(intStartBefore(i):intEndBefore(i));
            thisIntervalAfter = thisChannel(intStartAfter(i):intEndAfter(i));
            
            if length(thisIntervalBefore) == pointsPerWindow
                usableDataBefore(:, i) = thisIntervalBefore;
            else
                fprintf('Unexpected length for before-window at %d: Found %d, expected %d\n', ...
                    i, length(thisIntervalBefore), pointsPerWindow);
            end
            
            if length(thisIntervalAfter) == pointsPerWindow
                usableDataAfter(:, i) = thisIntervalAfter;
            else
                fprintf('Unexpected length for after-window at %d: Found %d, expected %d\n', ...
                    i, length(thisIntervalAfter), pointsPerWindow);
            end
        end
        
        fprintf('Size of usableDataBefore: %d %d\n', size(usableDataBefore,1), size(usableDataBefore,2));
        fprintf('Size of thisData: %d %d %d\n', size(thisData,1), size(thisData,2), size(thisData,3));

        thisData(:, activeIdx, 1:2:end) = usableDataBefore;
        thisData(:, activeIdx, 2:2:end) = usableDataAfter;
        
    end

    % Update the labels
    totalWindows = nWindowsParsed + nWindows;
    fileIdx = nWindowsParsed+1:totalWindows;
    labels.allWindows.mouse(fileIdx) = repmat({mousename}, nWindows, 1);
    labels.allWindows.expDate(fileIdx) = repmat({date}, nWindows, 1);
    labels.allWindows.time(fileIdx) = 1:nWindows;

    nWindowsParsed = totalWindows;
    dataCells = cat(3, dataCells, {thisData});

    fprintf('Day %s of %s loaded. %.1f minutes elapsed\n', date, mousename, toc/60);
end

% Concatenate all files
data = cell2mat(dataCells);

% Fill in remaining labels and save
labels.fsRaw = fs;
labels.windowLength = windowLength;

% Save data
save(saveFile, 'data', 'labels', '-v7.3');
end
