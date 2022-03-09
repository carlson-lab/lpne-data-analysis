datafile = 'cSDS_EPM.mat';
load(datafile, 'labels')

recordingOrder = zeros(size(labels.allWindows.expDate));

% get list of unique mouse names
mouseLabels = labels.allWindows.mouse; %mouse name for each window
uniqueMouseLabels = unique(mouseLabels); %list of unique mouse names

for mousename = uniqueMouseLabels %for each mouse
    
    mouseIdx = strcmp(mouseLabels, mousename); %the indices of windows that the current mouse has data in
    
    % get dates
    dateLabels = labels.allWindows.expDate(mouseIdx);
    date = unique(dateLabels);
    disp(date)

    % assuming 2 dates per mouse, one pre-manipulation, one post-manipulation
    %also assuming that the earlier date is 
    %first_date = date(1)



    if length(date) == 1
        date1 = date(1);
        date2 = 'NODATE';

    elseif length(date) == 2
        date1 = date(1);
        date2 = date(2);
    end
    



    

    % create new label subbing dates out for date index
    
    for d = 1:length(date)
        dateIdx = strcmp(labels.allWindows.expDate, date(d));
        recordingIdx = mouseIdx & dateIdx;
        
        if strcmp(date{d}, date2)
            recordingOrder(recordingIdx) = 2;
        elseif strcmp(date{d}, date1)
            recordingOrder(recordingIdx) = 1;
        else
            error('unrecognized date!')
        end
    end
end

labels.allWindows.recordingOrder = recordingOrder;

group = cell(1,length(mouseLabels));

for i=1:length(mouseLabels)
    group{i} = ['Delirium'];
end
labels.allWindows.groups = group;



save(datafile, 'labels', '-append')

