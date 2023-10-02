# lpne-data-analysis, frame_windows branch

The frame windows branch allows one to slice the LFPs into windows to the sides of specified timepoint. You can choose the amount of windows that you want to create to the either side of each timepoint. The way to input the timepoints for slicing is to make .mat files in the standard naming scheme for each LFP file, but instead of the file ending in _LFP.mat, it will end in _FRAME.mat. All of these files should be in a folder called 'Frame', which itself is contained in the same folder that also contains the Data and CHANS folders as well. The contents of this file should be a #_timepoints by 1 double that contain all of the timepoints. The variable needs to be named 'frames' in order to be read in. This branch automatically formats windows in this manner; you do not need to specify that you want to do this. You need to use another branch if you don't want to slice windows in this way. Preprocessing the windows and making the features are unchanged.



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
