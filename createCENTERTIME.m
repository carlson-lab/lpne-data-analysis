prompt = 'How long is the largest window of features you are going to create (in seconds)?' ;
largest_window_size = input(prompt);


prompt = 'How far apart do you want each timepoint to be (in seconds)? ';
timepointSepLen = input(prompt);

prompt = 'What is the sampling frequency?';
sampleFreq = input(prompt);

prompt = 'true or false: you have recordings stitched together that are all the same length?';
stitched = input(prompt);

if stitched
    prompt = 'How long are the recordings of each file you stitched together?';
    recLen = input(prompt);

    %read in LFP file
    projectFolder = uigetdir('.', 'Select folder containing Data & CHANS subfolders ');
    dataFolder = [projectFolder '/Data/'];
    centerFolder = [projectFolder '/CENTER_TIME/'];
    
    dataList = dir([dataFolder '*_LFP.mat']);
    
    nSessions = length(dataList);
    for k = 1:nSessions
        thisFile = dataList(k);
        if thisFile.isdir, continue, end
    
        % clear channel data from last file and load it from the next file
        clear('-regexp','_\d\d')
        filename = thisFile.name;
        dataPath = [dataFolder filename];
        dataFile = matfile(dataPath);
        varsList = who(dataFile);
        variableName = char(varsList(1));
        variable = load(dataPath, variableName);
        varSize = size(variable.(variableName));
        lenLFP = varSize(2);
    
        windowSample = largest_window_size * sampleFreq;
        startPoint = windowSample/2;
    
        sepSample = timepointSepLen * sampleFreq;
    
        timepoints = [];
        currentPoint = startPoint;
        recSample = recLen*sampleFreq;
        

        while currentPoint < lenLFP - startPoint
            if rem(currentPoint, recSample) ~= 0
                timepoints = [timepoints ; currentPoint];
                currentPoint = currentPoint + sepSample;
            else
                disp(currentPoint)
                currentPoint = currentPoint + sepSample;
            end
        end
    
       T = timepoints;
       percProgTraveledPath = zeros(length(T),1);
       trial = zeros(length(T),1);
    
       newFilename = strrep(filename, 'LFP', 'CENTER');
       save(strcat(centerFolder,newFilename), 'T', 'percProgTraveledPath', 'trial')
       
    end
    disp('Done')
    clear all

else

    %read in LFP file
    projectFolder = uigetdir('.', 'Select folder containing Data & CHANS subfolders ');
    dataFolder = [projectFolder '/Data/'];
    centerFolder = [projectFolder '/CENTER_TIME/'];
    
    dataList = dir([dataFolder '*_LFP.mat']);
    
    nSessions = length(dataList);
    for k = 1:nSessions
        thisFile = dataList(k);
        if thisFile.isdir, continue, end
    
        % clear channel data from last file and load it from the next file
        clear('-regexp','_\d\d')
        filename = thisFile.name;
        dataPath = [dataFolder filename];
        dataFile = matfile(dataPath);
        varsList = who(dataFile);
        variableName = char(varsList(1));
        variable = load(dataPath, variableName);
        varSize = size(variable.(variableName));
        lenLFP = varSize(2);
    
        windowSample = largest_window_size * sampleFreq;
        startPoint = windowSample/2;
    
        sepSample = timepointSepLen * sampleFreq;
    
        timepoints = [];
        currentPoint = startPoint;
        while currentPoint < lenLFP - startPoint
            timepoints = [timepoints ; currentPoint];
            currentPoint = currentPoint + sepSample;
        end
    
       T = timepoints;
       percProgTraveledPath = zeros(length(T),1);
       trial = zeros(length(T),1);
    
       newFilename = strrep(filename, 'LFP', 'CENTER');
       save(strcat(centerFolder,newFilename), 'T', 'percProgTraveledPath', 'trial')

       
    end
    disp('Done')
    clear all
end