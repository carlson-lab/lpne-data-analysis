function [averagedData, windowTimes ] = removeNaNs(averagedData, windowTimes, deleteWin, S)
% Removes any remaining NaN values from averagedData, creating a new
% segment for any discontinuous windows.

% First finds the cells to remove in each segment. Then, iterates by
% segment. For each time through the loop:
% Notes the position of the current cell and the next cell. 
% Pushes everything after the current cell backwards by the number of
% deleted windows / data to be added
% Updates the new next cell position (ie if it was 2, then 5 new cells were
% added, it would now be 7)
% Working backwards through available empty windows and through windows to
% delete, copies everything after the last remaining deletable window to
% the last available cell
% Removes everything including and following that last deletable window
% from the starting cell
% Does this process again for timeWindows - move everything after the
% delete to the last available, deletes everything including and after the
% delete from the original.
% Move to the next original cell and repeat
%
% After this has finished for each cell, find all empty cells(ie if you
% were deleting two windows adjacent, there would be nothing to move
% between them and so the cell would be empty), and delete them from data
% and timewindows

% INPUT:
% averagedData: Cell array of segments, each one of which is N by M where N is 
%   the number of time points, and M is the number of channels. Contiguous data is 
%   in a single segment
% windowTimes: Cell array of start and end times for each window within a channel. Each 
%   cell array is 2 by W, where W is the number of windows per channel.
% deleteWin: cell array containing indices of windowTimes that contain NaN
%   and therefore the windows need to be removed
% S: Number of original segments
%
% OUTPUT:
% averagedData: Same as before, but without any NaNs, and with each
%   discontinuous set of windows in its own segment
% windowTimes: same as before, but matched to the edited averagedData

% Step 1: Find the number of new cells to make
newSegsWin = {};

for s=1:S
    deleteWin{s} = unique(deleteWin{s});
    newSegsWin{s} = length(deleteWin{s});
end

% Initialize variables to help with indexing
startCell = 1;
nextStart = 2;

for s = 1:S
    % Step 2: Add all empty cells in the right places, after the cell
    % containing their data
    
    for j = length(averagedData):-1:nextStart
        averagedData{j+newSegsWin{s}} = averagedData{j};
        windowTimes{j+newSegsWin{s}} = windowTimes{j};
    end
    
    deleteCell = unique(deleteWin{s});
    
    % track the position of the next unaltered cell
    nextStart = nextStart + newSegsWin{s};
    
    % Step 3: Copy data into new cells, backwards, and remove duplicate data from original
    % cell
    fillPos = nextStart;
    for delIndex = length(deleteCell):-1:1
        % Copy only after the nan index, not inclusive
        
        % if the deleted index isn't the last window in the segment,
        if deleteCell(delIndex)<length(windowTimes{startCell})
            % Working backward, copy the chunk of data (not including the NaN) to the farthest away
            % open cell. Use the index to access the start point from
            % windowTimes, being sure to normalize if windowTimes doesn't
            % start at 1.
        averagedData{fillPos - 1} = averagedData{startCell}(windowTimes{startCell}(1,deleteCell(delIndex)+1)-windowTimes{startCell}(1)+1:end,:);
        end
        % remove from the original cell everything after the nan index,
        % including that window
        averagedData{startCell}(windowTimes{startCell}(1,deleteCell(delIndex))-windowTimes{startCell}(1)+1:end,:) = [];
    %the next chunk of data is closer to the starting cell
        fillPos = fillPos - 1;
        
    end
    %reset fillPos to repeat this loop for timeWindows. They're not done
    %simultaneously because the index from timeWindows is needed, and
    %therefore can't be altered until averagedData is complete.
    fillPos = nextStart;
    % repeat the loop for timeWindows
    for delIndex = length(deleteCell):-1:1
        if deleteCell(delIndex)<length(windowTimes{startCell})
        windowTimes{fillPos - 1} = windowTimes{startCell}(:, deleteCell(delIndex)+1:end);
        end
        windowTimes{startCell}(:,deleteCell(delIndex):end) = [];
        % move back a cell for filling
        fillPos = fillPos-1;
    end
    
    % Iterate forward to set the next start point,and next unaltered cell
    startCell = nextStart;
    nextStart = nextStart +1;
end

% Find the empty values, ie from where two adjacent windows contained NaNs
rmWin = ~cellfun('isempty',averagedData);
% remove them from both windowTimes and averagedData.
windowTimes = windowTimes(rmWin);
averagedData = averagedData(rmWin);
    