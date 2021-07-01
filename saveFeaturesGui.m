function saveFeaturesGui(saveFile)
% estimate spectral features
%
% In previous versions, there was an additional input parameter called
% "options." In this version a GUI pops up and asks for these options from
% the user.
%
% INPUTS
% saveFile: name of '.mat' file containing the data and labels
%   variables. And to which the processed data will be saved
%
% LOADED VARIABLES
% X: Preprocessed (filtered, averaged, checked for saturation) data. NxAxW array. A is
%       the # of areas. N=number of frequency points per
%       window. W=number of time windows.
% labels: Structure containing labeling infomation for data
%   FIELDS USED
%   fsRaw: sampling frequency of unprocessed data (Hz)
%   area: cell array of labels for each area corresponding to
%       the second dimension of xFft
%   fs: sampling frequency of processed data (Hz)
%   windows: structure containing relevant labels pertaining to
%           individual windows. Each field should be a vector / array with
%           one element corresponding to each window. Suggested fields:
%           date, etc. Must contain 'mouse', 'expDate', and 'time' fields.
%   windowLength: length of windows (s)
%
% SAVED VARIABLES
% power: MxNxP matrix of power values where M is frequency, N is brain area,
%     and P is time window
% coherence: MxNxPxQ array to store coherence variables where M is frequency,
%     N is time window, and P and Q are the two brain areas where coherence is
%     calculated
% granger: PxFxW array to store granger causality values. P
%     iterates over directed pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% instant: PxFxW array to store instantaneous causality values. P
%     iterates over undirected pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% directedSpectrum: PxFxW array to store 'full' model directed spectrum
%     features. P iterates over directed pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% pwDirectedSpectrum: PxFxW array to store pairwise directed spectrum
%     features. P iterates over directed pairs of regions, F iterates over
%     frequencies, W iterates over windows.
% fft: fourier transform of X
% labels: See above, with
%   ADDED FIELDS
%   f: integer frequency of processed data
%   powerFeatures: MxN matrix of string labels describing the
%       features represented in power. M is frequency, N is
%       brain area.
%   cohFeatures: MxPxQ array of string labels describing the
%       features represented in coherence. M is frequency, P
%       and Q are the two brain areas where coherence is calculated.
%   gcFeatures: PxF array of string labels describing the
%       features represented in granger. P iterates over
%       directed pairs of regions, F iterates over frequencies.
%   instFeatures: PxF array of string labels describing the
%       features represented in instArray. P iterates over
%       undirected pairs of regions, F iterates over frequencies.
%   dsFeatures: PxF array of string labels describing the
%       features represented in directedSpectrum and pwDirectedSpectrum
%       P iterates over directed pairs of regions, F iterates over
%       frequencies.

%% Get options from GUI
fprintf('Make sure you are using matlab version R2019a or later')
myGui = gui();
while isvalid(myGui)
    options = myGui.getOptions();
    pause(0.001);
end 
if isvalid(myGui)
   % one last check on the off chance that someone clicked an option and closed the app 
   % in 1ms
   options = myGui.getOptions();
end

saveFeatures(saveFile, options)

end
