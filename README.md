# lpne-data-analysis, centered_windows_functionality branch

If a 'true' boolean is passed to the input useCenteredWindows, then you will need a 'CENTER_TIME' folder in the project folder containing files with a list of timestamps to create windows around. The files are named the same as the LFP files, but end in CENTER.mat instead of LFP.mat. Each center file contains three variables:

1. T, which contains the timestamps
2. trial, what trial is each timestamp occuring in
3. percProgTraveledPath, the percent progress through the task the mouse is at that timestamp

If 2 and 3 do not mean anything to you, you can make dummy variables of them.

A special note for formatWindows is that it asks "Enter the largest window length used in your analysis (s)". This comes from that this branch can make multiple lengths of features (i.e. 1 second power features and 2 second coherence features). So if your features are 1 second, 2 seconds, and 2 seconds, you would input 2 to this. You will specify how long each feature is later.

When using saveFeatures, we will have a struct called 'options' contain information about what length each feature is. The field in options should be called featureSizes.

featureSizes: (Optional) Requires featureList to use. Cell array of doubles indicating the window size of each feature in featureList.
The doubles should be in the same order as the features specified in featureList. Only use if centered windows earlier. If you don't know what this means, don't include this.

featureSizes must be used if windows are centered







### A detailed description of how to use the full pipeline with control over individual steps is given in *NMFDemo.ipynb*

The main steps involved are:
1. **Format Data** \[MATLAB\]
2. **Preprocessing** \[MATLAB\]
3. **Calculate frequency-based features (e.g. power)** \[MATLAB\]
4. **Prepare data in Python** \[Python\]
5. **Train and evaluate NMF model** \[Python\]
6. **Save results** \[Python\]
7. **Interpret results** \[Python & MATLAB\]

### To run all the MATLAB steps together with the default options, run the following line:
```matlab
runFeaturePipeline('saveFile.mat')
```
replacing *saveFile.mat* with whatever filename you wish to save your preprocessed data and features in.

If you are using this data for 'backprojecting' into an already-trained model, instead run the following:
```matlab
runFeaturePipeline('saveFile.mat',[],[],'trainFile.mat')
```
replacing *saveFile.mat* with your save location as above and *trainFile.mat* with the file containing the preprocessed training data for the original model.

See the comments located in *runFeaturePipeline.m* for info on running the steps 1-3 of the pipeline with non-default options, or see *NMFDemo.ipynb* for more details on running each step individually.

### Additional useful codes included:

* NMFDemo.ipynb (Demo for proper use of data_tools.py and validation_tools.py, edit as necessary)

* datautils.segments2array.m (used in saveFeatures to convert from segments back to array. Can also be used by hand).

* datautils.array2segments.m (used to convert arrays to cell format, generally not needed but good to have)

* channelSpectrums.mat (generates power spectrums of individual channels saved in a preprocessed data file)

* mrl_phase_lag_estimates (estimates the phase offset between two signals in a specified frequency band)
