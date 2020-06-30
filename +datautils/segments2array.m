function [data] = segments2array(dataSegments, windowTimes)
% Convert data from being stored as contiguous segments to as 3-D array of
% windows
% INPUT
% dataSegments: cell array with contiguous segments of 2D data MxN, where M
%    is frequency and N is brain area; each cell is a window
% windowTimes: cell array with start / end indices for each window; each
%   cell corresponds to a segment
%
% OUTPUT
% data: matrix MxNxP where M is frequency, N is brain area, and P is time
%   window

% accumulate windows from all segments into array
data = [];
for s = 1:numel(dataSegments)
  
  theseTimes = windowTimes{s};
  theseTimes = theseTimes - theseTimes(1)+1;
  W = size(theseTimes, 2);

  % extract data corresponding to each time window and concatenate to data
  % along 3rd dimension
  for w = 1:W
    thisWindow = dataSegments{s}(theseTimes(1,w):theseTimes(2,w),:);
    data = cat(3, data, thisWindow);
  end
  
end
end

