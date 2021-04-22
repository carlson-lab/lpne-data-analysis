function savePreprocessedDataSubset(loadFile, saveFile, subIdx)
origData = load(loadFile);

if numel(subIdx) == 1
    W = length(origData.labels.windows.expDate);
    nSamples = subIdx;
    subIdx = false(W,1);
    subIdx(randsample(W, nSamples)) = true;
end

newData.X = origData.X(:,:,subIdx);

if isfield(origData, 'power')
    newData.power = origData.power(:,:,subIdx);
    labels.powerFeatures = origData.labels.powerFeatures;
end
if isfield(origData, 'coherence')
    newData.coherence = origData.coherence(:,subIdx,:,:);
    labels.cohFeatures = origData.labels.cohFeatures;
end
if isfield(origData, 'granger')
    newData.granger = origData.granger(:,:,subIdx);
    labels.gcFeatures = origData.labels.gcFeatures;
end
if isfield(origData, 'directionality')
    newData.directionality = origData.directionality(:,:,subIdx);
    labels.ldFeatures = origData.labels.ldFeatures;
end
if isfield(origData, 'directionality_pairwise')
    newData.directionality_pairwise = origData.directionality_pairwise(:,:,subIdx);
    labels.ldFeatures = origData.labels.ldFeatures;
end
if isfield(origData, 'directionality_cond')
    newData.directionality_cond = origData.directionality_cond(:,:,subIdx);
    labels.ldcFeatures = origData.labels.ldcFeatures;
end
if isfield(origData, 'xFft')
    newData.xFft = origData.xFft(:,:,subIdx);
end

% subset labels
labels.channel = origData.labels.channel;
labels.channelArea = origData.labels.channelArea;
labels.fsRaw = origData.labels.fsRaw;
labels.windowLength = origData.labels.windowLength;
labels.area = origData.labels.area;
labels.s = origData.labels.s;
labels.fs = origData.labels.fs;
labels.f = origData.labels.f;

labels.windows = origData.labels.windows;
for f = fieldnames(labels.windows)'
    labels.windows.(f{1}) = labels.windows.(f{1})(subIdx);
end
newData.labels = labels;

save(saveFile, '-struct','newData', '-v7.3')