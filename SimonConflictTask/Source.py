import pandas as pd
import scipy.io as sio
import os
import numpy as np
import mne
import matplotlib.pyplot as plt

# Import datasets, classifiers and performance metrics
from sklearn import datasets, svm, metrics, linear_model, model_selection

# Load the .mat file
matlabFile = sio.loadmat("SimonConflictTask/Subject_ID_25/Data and design matrices/ID_25_epochsAll_designMat.mat")

# Turn the matlabFile into a Panda dataframe for easier manipulation
labels = pd.DataFrame(matlabFile["epochs_all"], columns=["epochN", "artifactFlag", "binIndicator", "avg_reward", "reward", "iti", 
                                                         "trial_num", "run_num", "trial_in_run", "rt", "congruencyDC", "accDC", 
                                                         "previousCongruencyDC", "previousAcc", "trialType", "keyRep", "cnvAmplitude"])

# Only consider trials with valid data by removing NaN and artifactFlag = 0 values
labels = labels[labels.artifactFlag != 0].dropna()

# Label we are trying to predict
y = labels.congruencyDC.real

# Load selected channels around the FCz area from the EEG data
# After loading the EEG data the 'eegData' variable is composed by trials[868],channels[9],time[2560]
eegData = mne.read_epochs_eeglab("SimonConflictTask/Subject_ID_25/Data and design matrices/ID_25_epochsAll.set").pick_channels(["F1","Fz","F2","FC1","FCz","FC2","C1","Cz","C2"]).get_data()

# Now we remove the trials from the eegData that have invalid data (as established before)
X = np.ndarray(shape=(labels.epochN.size, eegData.shape[1], eegData.shape[2]), dtype='float64')
counter = 0
for index in labels.epochN.astype(int):
    X[counter] += eegData[index - 1,:X.shape[1],:X.shape[2]]
    counter += 1

# keep track of the accuracy over time
predictorTime = pd.DataFrame(columns=["Time", "Accuracy"])

# We use a Support Vector Machine with a linear kernel
classifier = svm.SVC()

# Train our classifier for 32 timepoint in the EEGData
timePoints = 32
for time in range(0, X.shape[2], int(2560 / timePoints)):

    # Create a training and testing set of labels/X where 10% of the data is used for training
    X_train, X_test, y_train, y_test = model_selection.train_test_split(X[:X.shape[0],:X.shape[1],time], y, train_size=0.1)
        
    # Train the classifier
    classifier.fit(X_train, y_train)

    # Tries to classifies the remaing 90% of the breast cancer data set
    y_pred = classifier.predict(X_test)

    # keep track of our prediction over time
    predictorTime = predictorTime.append({'Time': time, 'Accuracy': metrics.accuracy_score(y_test, y_pred)}, ignore_index=True)

# Visualize the results
xData = predictorTime.Time.real
yData = predictorTime.Accuracy.real

# Plot the accuracy,time graph
plt.plot(xData, yData, "r.-")

# Plot the trendline
z = np.polyfit(xData, yData, 1)
p = np.poly1d(z)
plt.plot(xData,p(xData),"b--")

# Set the grpah info
plt.xlabel("Time")
plt.ylabel("Accuracy")
plt.title("Mean accuracy: {:.2%}".format(predictorTime.Accuracy.mean(axis=0)))
plt.show()