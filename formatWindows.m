function formatWindows(saveFile, projectFolder, chanFile, fs, windowLength)

% In addition to the standard Data and CHANS folders, we also need...
%... a Frame folder, which contains files ending in _FRAME.mat, with...
%... the same naming scheme as the Data and CHANS files.
% The FRAME files should contain an nx1 double of the frames you..
%... want to use to create the windows to the either side of.


% Fill in default arguments and get input from the user.
if nargin < 2
    projectFolder = uigetdir('.', 'Select folder containing Data & CHANS subfolders');
    dummy = input(['Make sure areas in channel info file match other ' ...
        'datasets you plan to combine with this one!!\n', 'ENTER to continue']);
end
if nargin < 3
    [chanFile, chanPath] = uigetfile([projectFolder '/*.xls*'], 'Select channel info file');
    chanFile = [chanPath chanFile];
end
if nargin < 4
    inputs = inputdlg({'Enter sampling rate (Hz):'});
    fs = str2double(inputs{1}); % sampling rate, Hz
end
if nargin < 5
    inputs = inputdlg({'Enter window length (s)'});
    windowLength = str2double(inputs{1}); % length of one window (s)
end
if nargin < 6
    inputs = inputdlg({'Enter the number of windows to create before each frame:'});
    numWindowsBefore = str2double(inputs{1});
end

if nargin < 7
    inputs = inputdlg({'Enter the number of windows to create after each frame:'});
    numWindowsAfter = str2double(inputs{1});
end

pointsPerWindow = fs*windowLength;

% Load channel info (and strip of unwanted ' symbols)
chanData = readtable(chanFile, 'ReadVariableNames', false);
chanData = strrep(chanData{:,:}, "'", '');
channames = chanData(:,1);
chanareas = chanData(:,2);
labels.channel = channames;
labels.channelArea = chanareas;

% Initialize labels as an empty struct if it does not already exist
if ~isfield(labels, 'allWindows')
    labels.allWindows = struct();
end

% Initialize labels.allWindows.frame as an empty array if it does not already exist
if ~isfield(labels.allWindows, 'frame')
    labels.allWindows.frame = [];
    labels.allWindows.mouse = [];
    labels.allWindows.expDate = [];
    labels.allWindows.time = [];
end


% For each recording file, load and slice data
chansFolder = [projectFolder '/CHANS/'];
dataFolder = [projectFolder '/Data/'];
dataList = dir([dataFolder '*_LFP.mat']);
nSessions = length(dataList);

% Initialize variables to be used in loops below
dataCells = {};
tic
nWindowsParsed = 0;

nWindows = numWindowsBefore + numWindowsAfter; %this line was added and is not in the original

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

    % Load the frame file for the current LFP file
    frameFile = [projectFolder '/Frame/' regexprep(thisFile.name, 'LFP.mat', 'FRAME.mat')];
    if exist(frameFile, 'file')
        load(frameFile, 'frames')
    else
        error('Frame file not found for %s', thisFile.name);
    end

    % Extract mouse name and experiment date from filename
    nameParts = split(thisFile.name,'_');
    mousename = nameParts{1};
    date = nameParts{2};
    timeVar = 0;
    % Create data structure for this file's windows
    thisData = NaN(pointsPerWindow, length(channames), nWindows * length(frames));
    disp(size(thisData))
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
        
        totalWindows = numWindowsBefore + numWindowsAfter;
        
        for f = 1:length(frames)
            frameTime = frames(f);
            
            windowCounter = 0;
            for i = 1:numWindowsBefore
                intStart = frameTime - (numWindowsBefore - i + 1) * windowLength * fs;
                intEnd = frameTime - (numWindowsBefore - i) * windowLength * fs - 1;
    
                thisInterval = [];  % Initialize to an empty array
                
                if intStart > 0 && intEnd <= length(thisChannel)
                    fprintf('Doing before window for frame %d: intStart=%d, intEnd=%d\n', f, intStart, intEnd);
                    thisInterval = thisChannel(intStart:intEnd);
                    windowCounter = windowCounter + 1;
                    if c == 1
                        labels.allWindows.frame = [labels.allWindows.frame; frameTime];
                        labels.allWindows.mouse = [labels.allWindows.mouse; mousename];
                        labels.allWindows.expDate = [labels.allWindows.expDate; date];
                        timeVar = timeVar + 1; 
                        labels.allWindows.time = [labels.allWindows.time; timeVar];

                    end
                else
                    % Either skip this window or print a warning message.
                    fprintf('Skipping before window for frame %d: intStart=%d, intEnd=%d\n', f, intStart, intEnd);
                end

                if length(thisInterval) == pointsPerWindow
                    thisData(:, activeIdx, (f-1)*totalWindows + windowCounter) = thisInterval;
                end
            end
            
            for j = 1:numWindowsAfter
                intStart = frameTime + (j - 1) * windowLength * fs + 1;
                intEnd = frameTime + j * windowLength * fs;
    
                thisInterval = [];  % Initialize to an empty array
                
                if intStart > 0 && intEnd <= length(thisChannel)
                    fprintf('Doing after window for frame %d: intStart=%d, intEnd=%d\n', f, intStart, intEnd);
                    thisInterval = thisChannel(intStart:intEnd);
                    windowCounter = windowCounter + 1;
                    if c == 1 
                        labels.allWindows.frame = [labels.allWindows.frame; frameTime];
                        labels.allWindows.mouse = [labels.allWindows.mouse; cellstr(mousename)];
                        labels.allWindows.expDate = [labels.allWindows.expDate; cellstr(date)];
                        timeVar = timeVar + 1; 
                        labels.allWindows.time = [labels.allWindows.time; timeVar];
                    end
                else
                    % Either skip this window or print a warning message.
                    fprintf('Skipping after window for frame %d: intStart=%d, intEnd=%d\n', f, intStart, intEnd);
                end
                

                if length(thisInterval) == pointsPerWindow
                    thisData(:, activeIdx, (f-1)*totalWindows + windowCounter) = thisInterval;
                end
            end
        end
    end



    % Find the indices of zeros in labels.allWindows.frame
    zeroIndices = find(labels.allWindows.frame == 0);

    % Remove those indices
    labels.allWindows.frame(zeroIndices) = [];

    % Find slices (along 3rd dimension) where all values are NaN
    allNanSlices = all(all(isnan(thisData), 1), 2);
    
    % Remove those slices
    thisData(:, :, allNanSlices) = [];
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
