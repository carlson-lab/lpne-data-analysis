function addDayLabel(saveFile)
% add labels giving the ordering of each day for recordings for each mouse
% assumes that the file has a labels structure and that
% labels.allWindows.expDate exists. labels.allWindows.expDate must have
% dates saved in 'mmddyy' format

load(saveFile,'labels')

mouseNames = unique(labels.allWindows.mouse); % get unique mouse names

labels.allWindows.dayOrder = zeros(size(labels.allWindows.expDate));
if isfield(labels, 'windows')
    labels.windows.dayOrder = zeros(size(labels.windows.expDate));
end

for m = 1:numel(mouseNames)
    labels.allWindows = addMouseDays(labels.allWindows, mouseNames{m});
    if isfield(labels, 'windows')
       labels.windows = addMouseDays(labels.windows, mouseNames{m}); 
    end
end

save(saveFile, 'labels', '-append')

end

function windows = addMouseDays(windows, mouseName)
    thisMouseIdx = find(strcmp(windows.mouse,mouseName));
    thisMouseDates = windows.expDate(thisMouseIdx);
    redate = cellfun(@(x) [x(5:6) x(1:4)],thisMouseDates,...
        'UniformOutput',false);
    day = sort(unique(redate));
    for d = 1:length(day)
        thisDay = strcmp(redate, day(d));
        windows.dayOrder( thisMouseIdx(thisDay)) = d;
    end
end