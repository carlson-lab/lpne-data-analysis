function [ dataSegments, timeWindows ] = array2segments(data, startTimes)
% Convert data from being stored  as 3-D array of
% windows to as cell array of contiguous segments. If optional input
% timeWindows is used, will combine contiguous windows together into a single segment.
%
% INPUT
% data: MxNxP matrix where M corresponds to the number of time points per
%     window, N is the number of channels, and P is the number of windows.
% startTimes: (Optional) vector with the start times of each window. If
%     used, will help combine contiguous windows into a single segment
% 
% OUTPUT
% dataSegments: cell list of QxR matrices where Q is the number of time
%     points and R is the number of channels. If timeWindows input is not
%     given, default is one window per segment. However, may contain >1
%     window if multiple windows are either overlapping or adjacent.
% timeWindows: cell list of start and end times for each segment relative to
%     its start time (ie [1 301; 300 600]).


% initialize variables
W = size(data, 3);
T = size(data, 1);
dataSegments = {};
timeWindows = {};

% if timeWindows is not given, no data is continuous. Initialize 
% tstartTimes if needed. 
if nargin<2
    startTimes = cell(W,1);
end
 
% Set start time to -T so it is guaranteed to not be 
% continuous with cell 0.
currentTime = -T;

% for each window, check if the difference between times is greater than
% one window length. If so, assign to a new window. If not, combine to old window.


for w = 1:W 

   oldTime = currentTime;
  
   % if timeWindows is not given as a variable, say that the windows are 
   % not continuous. Also, if no more start times are given, say that windows are
   % not continuous
   if isempty(currentTime) || length(startTimes)<w
       timeGap = T+1;
   % otherwise check time gap
   else
       currentTime = startTimes(w);
       timeGap = currentTime - oldTime;
   end
   % if the time gap is less than or equal to a full window, add on the
   % missing portion to the previous window
   if timeGap < T+1
       startPt = oldTime+T+1-currentTime;
       dataSegments{end} = cat(1,dataSegments{end},data(startPt:end,:,w));
       % update time window to go from 1 to two full windows minus overlap
       timeWindows{end} = cat(2, timeWindows{end},[startPt; size(dataSegments{end},2)]);
     
   % if the time gap is greater than a full window, create a new segment
   else
        dataSegments = [dataSegments data(:,:,w)];
        timeWindows = [timeWindows [1; T]];
   end
   
   
end

end

