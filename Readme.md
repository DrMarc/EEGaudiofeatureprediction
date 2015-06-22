# EEG audio feature extraction

## Overview
These Matlab scripts compute several audio features and saliency models on binaural recordings.

## Dependencies

### Audio features

- Two!Ears fronted
- Auditory Modelling toolbox
- Large time-frequency analysis toolbox

### Saliency

- Kayser et al. 2005
- Kaya and Elhilali 2014

## Data format

You need one or more files of EEG recording in BDF format (other formats are easily incorporated by replacing the calls to eeglab's *pop_biosig*). You also need binaural recordings of stimuli heard during the EEG recordings in broadcast WAV format (this is Microsoft WAV plus a chunk that can hold marker timestamps). Many handheld sound recorders save data in this format. Markers in the BDF and WAV files have to correspond (this can be achieved by simultaneously pressing marker buttons on the EEG amplifier and the audio recorder).
The file names should be the same, except for the extension; for instance *subject1_session1.bdf* and a corresponding *subject1_session1.wav*.

To extract the features (this may take several hours!):

	features = combineEEGaudiofeatures(EEGfile,audiofile);
	save('run1_features.mat','features');

*features* is a nested Matlab struct with fields that hold the extracted audio and saliency features and a copies of these feature vectors resampled to, and aligned with, the EEG recording.

The EEG data is processed by repeated calls to *EEGanalysis.m*:

	% 1. step: preprocessing and manual artifact rejection
	EEGanalysis('run1.bdf');
	% 2. step: ICA of pruned data (non-interactive)
	EEGanalysis('run1.bdf');
	% 3. step: Prune ICA activations
	EEGanalysis('run1.bdf');
	% 4. step: Merge data sets from same session by
	% selecting several bdf files that will be concatenated
	EEGanalysis('run1.bdf');
	% 5. step: ICA on merged data (non-interactive)
	EEGanalysis('run1.bdf');
	% 6. step: Select ICs
	EEGanalysis('run1.bdf');
	% 7. step: Save IC spectrogram features
	EEGanalysis('run1.bdf'); % NOT IMPLEMENTED YET

Each step saves a result eeglab *.set* file and text or mat files with misc information (such as a list of rejected frames, or selected ICs). *EEGanalysis.m* detects these files and automatically starts the next step. (Delete unwanted results files to redo a step.)

Then call `preprocess_features('run1.bdf');` to remove rejected frames from the extracted features (so that they align with the processed EEG data), concatenate runs and z-score the feature values. The result should be aligned features (*run1_featurevec.mat*) and IC activations (*run1_icaact.mat*) ready for machine learning.

## TODO

* Step 7 of *EEGanalysis.mat* (compute spectrogram of IC activations)
* Reducing the spectral resolution of all 2D features to max 16 bins.
* The machine learning portion...