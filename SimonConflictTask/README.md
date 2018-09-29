## Simon Conflict Task

## EEG dataset information

* subjects: 22
* all raw data are available [here](https://www.dropbox.com/sh/q19b0zr9fx92eaa/AABEYrHMimw48SPJveipS3B2a?dl=0)
* 70 channels/sensors
* 868 trials/epochs/repetitions in the original data
* task: Simon conflict task
  * See example Simon task (not the same as the one used in our EEG study but similar ehough!): https://www.labvanced.com/player.html?id=2865
* ignore these channels (either reference or eye channels): 'Iz' 'M1' 'M2' 'SO1' 'LO1' 'IO1' 'LO2'

## Design matrices

All subjects' design matrices are available [here](https://www.dropbox.com/sh/lcuchtk85z154pl/AABzHKYZLv6tKjM3RJb9UBKPa?dl=0). Use _epochsAll_designMat.mat, not _epochsClean_designMat.mat

When analysing data, include only these epochs/trials:

* artifactFlag == 1 (clean trials) # I reversed coded it.... 
* accDC == 0 (correct trials) # again I reversed coded it...
* trials without missing values in the design matrix

### Column information in design matrices

* epochN: epoch/trial/repetition number
* artifactFlag: whether epoch contains noise (0: noise [exclude these epochs], 1: no noise [**include** these epochs])
* avg_reward: computational modeling of reward shown in each epoch
* reward: reward shown on display on each epoch
* trial_num: trial number
* rt: reaction time
* congruencyDC: congruent (0) vs. incongruent trials (1)
* accDC: correct (0) vs. incorrect trials (1)
