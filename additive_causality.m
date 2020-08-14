function [acArray, acFeatures] = additive_causality(data, areaList, fs, options)
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

if nargin < 2
    % fill with default parameters
    options=[];
end
options=fillDefaultOpts(options);

data = permute( data, [2,1,3]);

% determine which frequencies to evaluate
F = floor(fs/2) + 1;
f = sfreqs(F-1, fs);

if options.separateWindows
    % loop over each window seperately
    [C,~,W] = size(data);
else
    % run through loop once using all windows
    W = 1;
    C = size(data,1);
    thisData = data;
end

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

% estimate model order (using a subset of the data to conserve memory)
if W > options.ordSampleMaxW
    sampledData = data(:,:, randsample(W, options.ordSampleMaxW));
else
    sampledData = data;
end
[~,BIC] = tsdata_to_infocrit(sampledData, options.maxOrder, []);
[~,order] = min(BIC);

if options.parCores, pp = parpool([2 options.parCores]); end
a = tic;
for w = 1:W
    if options.separateWindows
        thisData = data(:, :, w);
        fprintf('Starting window %d: %2.1fs elapsed\n', w, toc(a))
    end

    parfor (cp = 1:nCP, options.parCores)
        c1 = channelPairs(1,cp); c2 = channelPairs(2,cp);
        
        % estimate spectral Granger causality in both directions
        % and instantaneous Granger causality
        [ac12, ac21] = calculateSpectralAC(thisData([c1 c2],:,:), 1, 2, order, F-1);
        
        % fill in output array
        acArray12(1,cp,:,w) = ac12; acArray21(1,cp,:,w) = ac21;
        
        if w==1
            % save labels corresponing to each element of feature arrays above
            acFeatures12(1, cp,:) = cellfun(@(x) [areaList{c1} '->' areaList{c2} ' ' x], ...
                fStrings, 'UniformOutput', false);
            acFeatures21(1, cp,:) = cellfun(@(x) [areaList{c2} '->' areaList{c1} ' ' x], ...
                fStrings, 'UniformOutput', false);
        end
    end
    
end
if options.parCores, delete(pp), end

%% combine values down here
acArray = reshape([acArray12; acArray21], 2*nCP, F, W);
acFeatures = reshape([acFeatures12; acFeatures21], 2*nCP, F);
end

function opts = fillDefaultOpts(opts)
if ~isfield(opts,'maxOrder'), opts.maxOrder = 50; end
if ~isfield(opts,'separateWindows'), opts.separateWindows = true; end
if ~isfield(opts,'parCores'), opts.parCores = 0; end
if ~isfield(opts,'ordSampleMaxW'), opts.ordSampleMaxW = 1000; end    
end

function [fxy, fyx] = calculateSpectralAC(U,x,y,p,fres,regmode)
%% calculateSpectralAC
%
% Calculate _unconditional_ additive frequency-domain MVGC (spectral multivariate
% Granger causality) from time series data by "traditional" method (as e.g.
% in GCCA toolbox)
%
% Adapted from the MVGC toolbox developed by Lionel Barnett and Anil K. Seth, 2012
%
%% Syntax
%
%     f = calculateSpectralAC(U,x,y,p,fres,regmode)
%
%% Arguments
%
% _input_
%
%     U          multi-trial time series data
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     p          model order (number of lags)
%     fres       frequency resolution
%     regmode    regression mode (default as in 'tsdata_to_var')
%
% _output_
%
%     fyx, fxy   unconditional additive spectral Granger causality
%
%% Description
%
% Returns the _unconditional_ additive frequency-domain (spectral) GC
% from the variable |Y| to the variable |X|, and vice versa, in the time 
% series data |U|, for model order |p|. Also calculated is the instantaneous
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. Call |freqs = <sfreqs.html sfreqs>(fres,fs)|, where
% |fs| is the sampling rate, to get a corresponding vector |freqs| of
% frequencies on |[0,fs/2]|. The regression mode is set by the |regmode|
% parameter, which may be |'LWR'| or |'OLS'| (see <tsdata_to_var.html
% |tsdata_to_var|> for details and defaults).
%
%%
if nargin < 6, regmode = []; end % ensure 'tsdata_to_var' default

n = size(U,1);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)

h = fres+1;

assert(isempty(z),'conditional spectral MVGC not available in GCCA mode');

xy = [x y];

U = U(xy,:,:); % extract variables, rearrange

nx = length(x);
n = length(xy);
x = 1:nx;
y = nx+1:n;

owstate = warn_supp;
[A,SIG] = tsdata_to_var(U,p,regmode);        % full regression
warn_test(owstate,  'in full regression - data non-stationary or colinear?');
if warn_if(isbad(A),'in full regression - regression failed'), return; end % show-stopper!

owstate = warn_supp;
[S,H] = var_to_cpsd(A,SIG,fres);             % spectrum & transfer function
warn_test(owstate,  'in spectral calculation');
if warn_if(isbad(S),'in spectral calculation - calculation failed'), return; end % show-stopper!

Hxy = H(x,y,:);
SIGyx = SIG(y,y)-(SIG(y,x)/SIG(x,x))*SIG(x,y); % partial covariance

Hyx = H(y,x,:);
SIGxy = SIG(x,x)-(SIG(x,y)/SIG(y,y))*SIG(y,x);

fyx = zeros(1,h); fxy = zeros(1,h);
for k = 1:h
    fyx(k) = real(Hxy(:,:,k)*SIGyx*Hxy(:,:, k)');
    fxy(k) = real(Hyx(:,:,k)*SIGxy*Hyx(:,:, k)');
end
end