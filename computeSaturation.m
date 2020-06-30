function [saturatedPoint, dataOpts] = computeSaturation(data, labels, dataOpts)
%ComputeSaturation - Calculate saturation points for LFP data
%   Generates a logical matrix indicating points at which the data for each 
%   channel seems to be saturated. The current version finds a channel to be 
%   saturated if data is above some threshold, which may be different for 
%   different channels.
%   INPUTS
%   data: cell vector, where each cell contains MxNXP array of the data for each
%       delay length. M is the # of time points. N is the # of channels. P
%       is the # of windows
%   labels: Structure containing labeling infomation for data
%       FIELDS (used here)
%       Channels: cell arrays of labels for channels
%   dataOpts: Structure containing data processing options
%       FIELDS (used here)
%       satThresh: can be a numerical vector or a string. If a vector, each
%           element gives the saturation detection threshold for the 
%           corresponding area in the data. If a string, uses a method for
%           automatic detection of thresholds: 'MAD' indicates to used
%           median absolute deviation method to detect outlier points; 
%           'SD' indicates to use standard deviation.
%   OUTPUTS
%   saturatedPoint: cell vector, where each cell contains locical MxNXP array
%       of the data for each delay length indicating whether a data point is
%       considered saturated (i.e. there is a signal artifact).  M is the # of
%       time points. N is the # of channels. P is the # of windows
%   thresholds: vector of upper threshold values for tagging saturation in each
%       area.

% initialize things
if isfield( labels, 'area')
  area = labels.area;
else
  area = unique( labels.channelArea );
end
nAreas = numel( area );
[nSamples, nChannels, nWindows] = size(data);
missing = isnan( data );

saturatedPoint = false(nSamples, nChannels, nWindows);
satThresh = zeros(1,nAreas);
% For each area determine outiers
for a = 1:nAreas
    if isnumeric(dataOpts.satThresh)
        thisThresh = dataOpts.satThresh(a);
    else
        thisThresh = dataOpts.satThresh;
    end
    channelsInArea = ismember(labels.channelArea, area{a});
    [thisSaturation, satThresh(a)] = markAreaSaturation(data(:, channelsInArea, :),...
                                                        thisThresh);
    
    % handle missing data
    missingHere = missing(:,channelsInArea,:);
    thisSaturation( missingHere ) = true;
    
    saturatedPoint(:, channelsInArea, :) = thisSaturation;
end

dataOpts.satThresh = satThresh;
end

function [saturatedPoints, upperThreshold] =  markAreaSaturation(areaData, thresh)
FULL_WINDOW_THRESHOLD = 0.5;
NO_RECORDING_SD = 1e-2;
THRESH_RATIO = 30;

[nSamples, nChannels, nWindows] = size(areaData);

envelopeMax = abs(hilbert(areaData));

if nargin < 2 || strcmp(thresh,'SD')
    areaSD = std(envelopeMax(:));
    upperThreshold = 3*areaSD;
elseif strcmp(thresh,'MAD')
    areaMD = median(abs(envelopeMax(:)),'omitnan');
    upperThreshold = 5*areaMD;
elseif isnumeric(thresh)
    upperThreshold = thresh;
else
    error('Unrecognized threshold method')
end
lowerThreshold = upperThreshold/THRESH_RATIO;

areaData = reshape( abs(areaData), nSamples, []);
saturatedPoints = false(size(areaData));

aboveUpper = areaData > upperThreshold;
belowLower = areaData < lowerThreshold;

% Calculate first difference in time dimension of indicator matrix for
% points above the upper threshold. The start of contiguous regions above
% threshold should be equal to 1 in the D array. Add a 1 to the end of D
% array so that end of each window is also marked for inspection below.
buff1 = ones(1,size(areaData,2));
buff0 = 0*buff1;
buffAU = cat(1,buff0,aboveUpper,buff1);
D = diff(buffAU);
[upperPoint, upperWindow] = find(D == 1);

% For each period where the signal is over the upper threshold, find the
% nearest points before and after for which the signal goes below the lower
% threshold. Set all points between those two as saturated.
thisWindow = 0;
prevSatEnd = 1;
for s = 1:numel(upperWindow)
  thisUpperPoint = upperPoint(s);
  lastWindow = thisWindow;
  thisWindow = upperWindow(s);
  
% Check if this point is in new window. If so, reset saved duration and
% end point. If not, skip to next point if this point is already known
% to be saturated.
  if thisWindow ~= lastWindow
    prevSatDuration = 0;
    prevSatEnd = 1;
  elseif thisUpperPoint < prevSatEnd
    continue
  end
  
% Check if this point is the end of a window. If so, set SatEnd and
% SatStart so that 'SatStart:SatEnd' is empty unless SatStart is changed
% to an earlier point. Else, set saturation block endpoints as described
% above.
  if thisUpperPoint == nSamples+1
    satStart = nSamples+1;
    satEnd = nSamples;
  else
    windowBL = belowLower(:, thisWindow);
    blBefore = windowBL(1:thisUpperPoint);
    blAfter = windowBL(thisUpperPoint:end);
    blBefore(1) = 1;
    blAfter(end) = 1;
    satStart = find(blBefore,1,'last');
    satEnd = find(blAfter,1,'first') + thisUpperPoint-1;
  end
  
% If the space between last and current saturation block is less than
% the average duration of the two blocks, consider space in between as
% saturated as well. If this is the case, set previous saturation block
% duration to maximum of previous or current.
  satDuration = satEnd - satStart + 1;
  if satStart-prevSatEnd < max(satDuration,prevSatDuration)
    satStart = prevSatEnd;
  end
  
  saturatedPoints(satStart:satEnd, thisWindow) = true;
  prevSatEnd = satEnd;
  prevSatDuration = satDuration;
end

% If channel is saturated for too much of a winddow, count whole window
% as saturated.
satFrac = sum(saturatedPoints, 1)/nSamples;
fullySatWindows = satFrac > FULL_WINDOW_THRESHOLD;
satWindowPts = repmat(fullySatWindows,nSamples,1,1);
saturatedPoints(satWindowPts) = true;

% mark windows where recording is stopped as saturated
noRecording = std(areaData,1) < NO_RECORDING_SD;
noRecording = repmat(noRecording, nSamples,1,1);
saturatedPoints(noRecording) = true;

saturatedPoints = reshape( saturatedPoints, nSamples, nChannels, ...
                           nWindows );
end
