function saveDataSubset(loadFile, saveFile, subIdx)
origData = load(loadFile);

if numel(subIdx) == 1
    W = length(origData.labels.windows.expDate);
    nSamples = subIdx;
    subIdx = false(W,1);
    subIdx(randsample(W, nSamples)) = true;
end

newData.data = origData.data(:,:,subIdx);

% subset labels
labels.channel = origData.labels.channel;
labels.channelArea = origData.labels.channelArea;
labels.fsRaw = origData.labels.fsRaw;
labels.windowLength = origData.labels.windowLength;

labels.allWindows = origData.labels.allWindows;
for f = fieldnames(labels.allWindows)'
    labels.allWindows.(f{1}) = labels.allWindows.(f{1})(subIdx);
end
newData.labels = labels;

save(saveFile, '-struct','newData', '-v7.3')