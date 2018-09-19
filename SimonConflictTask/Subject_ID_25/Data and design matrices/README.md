## Simon Conflict Task (Subject ID_25)

## EEG dataset information

* Subject ID: ID_25
* Raw data are available [here](https://www.dropbox.com/sh/g4rcelhy1mcok5t/AABhXFLr4hpQBTnYNO3XTLcAa?dl=0)
* 70 channels/sensors
* 868 trials/epochs/repetitions in the original data
* Task: Simon conflict task
  * See example Simon task (not the same as the one used in our EEG study but similar ehough!): https://www.labvanced.com/player.html?id=2865
* Basic results are contained in the Figures directory

EEG data single-trial information (information associated with each trial/epoch will be contained in these .mat or .txt files; for creating design matrices; we'll discuss these eventually...)

* ID_25_epochsAll_designMat.mat (in general, use this version)
* ID_25_epochsAll_designMat.txt (same as the .mat above)
* ID_25_epochsClean_designMat.mat (ignore this one)
* ID_25_epochsClean_designMat.txt (ignore this one)

## Useful/relevant variable/column information in the design matrices

* epochN: epoch/trial/repetition number

* artifactFlag: whether epoch contains noise (0: noise [exclude these epochs], 1: no noise [include these epochs])

* avg_reward: computational modeling of reward shown in each epoch

* reward: reward shown on display on each epoch

* trial_num: trial number

* rt: reaction time

* congruencyDC: congruent (0) vs. incongruent trials (1)

* accDC: correct (0) vs. incorrect trials (1)

## EEG time-domain dataset

MATLAB/EEGLAB; data are in the EEG.data structure: 70x2560x868 matrix: channel_time_trial) (MNE python should read these...

* ID_25_epochsAll.fdt
* ID_25_epochsAll.set

## EEG time-frequency dataset (MATLAB) (scipy.io should read .mat files)

* ID_25_singleTrials.mat
  * contains 4D matrix with fourier coefficients (complex numbers)
  * only 14-channels/sensors included
  * data are in the tf_fourierSpec matrix (547x14x40x103 matrix: trial_channel_frequency_time)