function [acArray, acFeatures, S] = additive_causality2(data, areaList, fs, f,...
    windowSize, parCores)
% Additive frequency-domain granger causality
%    Estimates additive frequency-domain granger causality from x->y and y->x. Estimates
%    are given for all frequency bands given by f.
%
%    INPUTS
%    X: Preprocessed (filtered, averaged, checked for saturation) data. NxAxW array. A is
%        the # of areas. N=number of frequency points per
%        window. W=number of time windows.
%    fs: Sampling rate (Hz)
%    areaList: ordered list of brain areas from which data was recorded
%    options: structure of options for this function. If any not included,
%        will be filled in with default values.
%    FIELDS
%      separateWindows: boolean indicating whether to model
%          windows seperately [Default: true]
%      maxOrder: maximum model order for model order estimation [Defaut: 50]
%      parCores: integer indicating number of cores to use for%
%          parallel computing. If 0, all tasks executed in serial. [Default: 0]
%      ordSampleMaxW: maximum number of windows to use in estimating the
%          order of the AR models used to calculate Granger values.
%    OUTPUTS
%    acArray: spectral granger causality estimates
%    acFeatures: PxF array of string labels describing the
%        features represented in acArray. P iterates over
%        directed pairs of regions, F iterates over frequencies.
%
% This function requires the MVGC toolbox. Portions of this code were taken
% from mvgc_demo script of that toolbox.

if nargin < 6
   parCores = 0; 
end

% determine which frequencies to evaluate
f = [0, f];
F = numel(f);

[~,C,W] = size(data);
channelPairs = combnk(1:C, 2)';

% number of directed and undirected pairs
nCP = size(channelPairs, 2);

% generate list of strings for frequencies
fStrings = compose('%d', f)';

% initialize arrays to be filled in loop below
acArray12 = zeros(1, nCP, F, W, 'single');
acArray21 = zeros(1, nCP, F, W, 'single');
acFeatures12 = cell(1, nCP, F);
acFeatures21 = cell(1, nCP, F);
if nargout > 2
    S = zeros(C,C,F-1,W);
end
    
if parCores, pp = parpool([2 parCores]); end
a = tic;
for w = 1:W
    thisData = data(:, :, w);
    fprintf('Starting window %d: %2.1fs elapsed\n', w, toc(a))
    
    thisCpsd = cpsd(thisData, thisData, windowSize, [], f, fs, 'mimo');
    thisCpsd = permute(thisCpsd, [2,3,1]);
    
    parfor (cp = 1:nCP, parCores)
    %for cp = 1:nCP  
        c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
        
        % estimate noise covariance and transfer matrix from pair cpsd
        pairCpsd = thisCpsd([c1,c2],[c1,c2],:);
        [ac12, ac21] = additive_causality_nonpar(pairCpsd, fs, F);
        
        % fill in output array
        acArray12(1,cp,:,w) = ac12; acArray21(1,cp,:,w) = ac21;
    end
    
    if nargout > 2
        % remove DC term from cpsd before saving
        S(:,:,:,w) = thisCpsd(:,:,2:end);
    end
    
    if w==1
        for cp = 1:nCP
            c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
            
            % save labels corresponing to each element of feature arrays above
            acFeatures12(1, cp,:) = cellfun(@(x) [areaList{c1} '->' areaList{c2} ' ' x], ...
                fStrings, 'UniformOutput', false);
            acFeatures21(1, cp,:) = cellfun(@(x) [areaList{c2} '->' areaList{c1} ' ' x], ...
                fStrings, 'UniformOutput', false);
        end
    end
    
end
if parCores, delete(pp), end

%% combine values down here
acArray = reshape([acArray12; acArray21], 2*nCP, F, W);
acFeatures = reshape([acFeatures12; acFeatures21], 2*nCP, F);
end

function [ac12, ac21] = additive_causality_nonpar(S, fs, h)
% calculate pairwise additive causality using non-parametric estimator

[H, SIG] = wilson_sf(S, fs);
        
H12 = H(1,2,:);
SIG21 = SIG(2,2)-(SIG(2,1)/SIG(1,1))*SIG(1,2); % partial covariance

H21 = H(2,1,:);
SIG12 = SIG(1,1)-(SIG(1,2)/SIG(2,2))*SIG(2,1);

ac21 = zeros(1,h); ac12 = zeros(1,h);
for k = 1:h
    ac21(k) = real(H12(:,:,k)*SIG21*H12(:,:, k)');
    ac12(k) = real(H21(:,:,k)*SIG12*H21(:,:, k)');
end
end