# lpne-data-analysis

Generic pre-processing and analysis code for LPNE data science. Takes in data from individual recordings, preprocesses and extracts features, then creates predictive models for tasks and tests performance. Main files are formatWindows.m, preprocessData.m, saveFeatures.m, data_tools.py and validation_tools.py.

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
