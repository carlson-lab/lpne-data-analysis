%function channelSpectrums(datafile)
% generates power spectrums of individual channels

%load(datafile,'xFft','labels','dataOpts')
s = labels.s;

% since lowFreq and highFreq are no longer options, shouldn't this be
% removed? and we should just use the whole data? 
freqUsed = s>dataOpts.lowFreq & s<dataOpts.highFreq;
xFft2 = xFft(freqUsed,:,:);
s2 = s(freqUsed);

% C is the number of channels 
C = size(xFft2,2);
fC = factor(C); % number of factors 
% if it's a prime number, factor C+1
if numel(fC) == 1, fC = factor(C+1); end
% a is the product of every 2 factors; b is the product of everything else
a = prod(fC(end:-2:1));
b = prod(fC(end-1:-2:1));

% for each channel
for c = 1:C
    % open that channel's subplot
    subplot(a,b,c)
    %open the fourier transform for all time for this channel
    thisChan = xFft2(:,c,:);
    % square these values, take the average according to windows, cut in
    % half - generating power spectral density
    psdx = mean(abs(thisChan).^2,3)/2;
    % plot a semilog plot with frequency space labels, then those times the
    % computed power spectral density - why not Welch?
    semilogx(s2,psdx.*s2')
    % label the y by this channel
    ylabel(labels.area{c},'Interpreter','none')
   % set x from start of frequency to end
    xlim([s2(1) s2(end)])
    yticks([])
    % set bounds
    xticks([2:2:10 20:20:100])
    xticklabels([2:2:10 20:20:100])
end
% give an overall title 
xlabel('Freq (Hz)')
suptitle('Power x Frequency')
