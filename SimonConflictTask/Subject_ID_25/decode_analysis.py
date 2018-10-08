import pandas as pd
import numpy as np
import os
import copy
import re
from scipy.io import loadmat
import matplotlib.pyplot as plt
import mne
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

pd.set_option('display.max_colwidth', 20, 'display.expand_frame_repr', False, 'display.max_rows', 16)
np.set_printoptions(edgeitems=10, infstr='inf', linewidth=150, nanstr='nan', precision=4, suppress=True, threshold=1000, formatter=None)


# set up directory
PATHS = {'CURRENT_DIR': os.getcwd()}
PATHS['DATA_DIR'] = os.path.join(PATHS['CURRENT_DIR'], 'Data and design matrices')
PATHS['TIME_DOMAIN_DATA_FILE'] = os.path.join(PATHS['DATA_DIR'], 'ID_25_epochsAll.set')
PATHS['DESIGN_MATRIX_FILE'] = os.path.join(PATHS['DATA_DIR'], 'ID_25_epochsAll_designMat.mat')


def load_design_matrix(design_mat_path, design_mat_name, header_name):
    """load design matrix and header, and convert to pandas dataframe"""
    mat = loadmat(design_mat_path)
    columns = np.hstack(mat[header_name].flat)
    df = pd.DataFrame(data=mat[design_mat_name], columns=columns)

    return df

df_design = load_design_matrix(PATHS['DESIGN_MATRIX_FILE'], 'epochs_all', 'epochs_allVars')

def cast_types(df, col_types):
    """convert column types"""
    for col in col_types:
        df[col] = df[col].astype(col_types[col])

    return df


df_design
col_types = {
    'epochN': float,
    'artifactFlag': 'category',
    'binIndicator': 'category',
    'avg_reward': float,
    'reward': float,
    'iti': float,
    'trial_num': float,
    'run_num': float,
    'trial_in_run': float,
    'rt': float,
    'congruencyDC': 'category',
    'accDC': 'category',
    'previousCongruencyDC': 'category',
    'previousAcc': 'category',
    'trialType': 'category',
    'keyRep': 'category',
    'cnvAmplitude': float
}

df_design = cast_types(df_design, col_types)
df_design.dtypes

epochs = mne.io.read_epochs_eeglab(PATHS['TIME_DOMAIN_DATA_FILE'], verbose=True)
epochs.metadata = df_design


def drop_na_epochs(epochs):
    """drop NaN values from epochs object and design matrix dataframe"""
    print("Original number of epochs in data and design matrix matches: {0} in data and {1} in design matrix".format(len(epochs), epochs.metadata.shape[0]))
    epochs.drop([any(x) for x in epochs.metadata.isnull().values])
    epochs.metadata = epochs.metadata.dropna()

    if len(epochs) != epochs.metadata.shape[0]:
        print("number of epochs in data and design matrix DO NOT MATCH: {0} in data and {1} in design matrix".format(len(epochs), epochs.metadata.shape[0]))
    else:
        print("New number of epochs in data and design matrix matches: {0}".format(len(epochs)))

    return epochs


drop_na_epochs(epochs)


epochs.filter(l_freq=None, h_freq=30, n_jobs=-1)  # filter (remove high-freq > 30 Hz)
epochs.crop(-0.2, 0.8)  # trim time
epochs.resample(sfreq=100, n_jobs=-1)  # resample to 100 Hz


def select_epochs(epochs, bools_to_keep):
    """select (drop) data/epochs based on boolean"""
    print("Original number of epochs in data and design matrix matches: {0} in data and {1} in design matrix".format(len(epochs), epochs.metadata.shape[0]))

    epochs.drop(np.logical_not(bools_to_keep))  # drop epochs

    if len(epochs) != epochs.metadata.shape[0]:
        print("number of epochs in data and design matrix DO NOT MATCH: {0} in data and {1} in design matrix".format(len(epochs), epochs.metadata.shape[0]))
    else:
        print("New number of epochs in data and design matrix matches: {0}".format(len(epochs)))

    return epochs


# select epochs and channels
epochs_subset = copy.deepcopy(epochs)  # deepcopy object so changes aren't made in place from now on
epochs_subset = select_epochs(epochs_subset, epochs_subset.metadata['artifactFlag'] == 1)  # only epochs without artifact
epochs_subset = select_epochs(epochs_subset, epochs_subset.metadata['accDC'] == 0)  # only accurate epochs
epochs_subset.drop_channels(['M1', 'M2', 'SO1', 'LO1', 'IO1', 'LO2'])
# epochs_subset.pick_channels(['FCz'])
epochs_subset.ch_names
epochs_subset.info


def set_events(epochs, target, event_ids):
    """update events and event_id"""
    epochs.events[:, 2] = epochs.metadata[target]
    epochs.event_id = event_ids
    print("Updated epochs.events, epochs.event_id")

    return epochs


def get_classifier_training_data(epochs, target, event_ids):
    """get features and targets for classification"""
    set_events(epochs, target, event_ids)
    epochs, _ = epochs.equalize_event_counts(sorted(event_ids))

    X = epochs.get_data()  # get data (features)
    y = epochs.events[:, 2]  # get target
    print("target categories: {0}".format(np.unique(y, return_counts=True)))

    return X, y


epochs_subset.metadata['congruencyDC'].value_counts()
X, y = get_classifier_training_data(epochs_subset, "congruencyDC", {'congruent': 0, 'incongruent': 1})
# X, y = get_classifier_training_data(epochs_subset, "accDC", {'correct': 0, 'incorrect': 1})


def decode(X, y, estimator, scoring, cv=10, n_jobs=-1):
    pipeline = make_pipeline(StandardScaler(), estimator)
    slider = mne.decoding.SlidingEstimator(pipeline, scoring=scoring)
    # slider.decod.fit(X, y)
    scores = cross_val_multiscore(slider, X, y, cv=cv, n_jobs=n_jobs)

    return slider, scores


def plot_scores(times, scores, scoring, title, chance=None, x_label='Time (s)'):
    """Plot averaged cross-validation scores"""

    if np.ndim(scores) > 1:
        scores = np.mean(scores, axis=0)  # mean cross-validation score

    fig, ax = plt.subplots()
    ax.plot(times, scores)
    ax.set_xlabel(x_label)
    ax.set_ylabel(scoring)

    ax.axvline(.0, color='k', linestyle='-')
    if chance is not None:
        ax.axhline(chance, color='k', linestyle='--', label='chance')
        ax.legend()

    ax.set_title(title)

    plt.show()


from sklearn import (linear_model, svm)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from mne.decoding import (SlidingEstimator, GeneralizingEstimator, cross_val_multiscore, LinearModel, get_coef)

decoder, scores = decode(X, y, LinearModel(LogisticRegression(C=1)), 'roc_auc', cv=10)
# decoder, scores = decode(X, y, svm.SVC(), 'roc_auc', cv=20)
scores.shape
scores
plot_scores(epochs_subset.times, scores, 'AUC', 'Decoding congruency (logistic regression)', chance=.5)


# TODO (1) turn into function; (2) save spatial filters/patterns; (3) use parallel loops
n_permutes=100
y_shuffled = copy.deepcopy(y)
scores_permuted = np.zeros((n_permutes, 10, 100))
scores_permuted.shape
p=0
for p in range(0, n_permutes):
    print("Running {0}th permutation...".format(p))
    np.random.shuffle(y_shuffled)  # shuffle target
    _, scores_temp = decode(X, y_shuffled, LinearModel(LogisticRegression(C=1)), 'roc_auc', cv=10)  # cross validate
    scores_permuted[p,...] = scores_temp  # save results
scores_permuted.shape
scores_permuted_mean = np.mean(scores_permuted, axis=1)
scores_permuted_mean = np.mean(scores_permuted_mean, axis=0)
scores_permuted_mean.shape
scores_permuted_sd = np.std(scores, axis=0)
scores_permuted_z = (np.mean(scores, axis=0) - scores_permuted_mean) / scores_permuted_sd

plot_scores(epochs_subset.times, scores_permuted[3,...], 'AUC', 'Decoding congruency (logistic regression)', chance=.5)
plot_scores(epochs_subset.times, scores_permuted_mean, 'AUC', 'Decoding congruency (logistic regression)', chance=.5)
plot_scores(epochs_subset.times, scores_permuted_z, 'Permuted z scores (AUC)', 'Decoding congruency (logistic regression)', chance=0.1)

# prepare a series of classifier applied at each time sample
clf = make_pipeline(StandardScaler(), LinearModel(LogisticRegression(C=1)))
# clf = make_pipeline(StandardScaler(), LinearModel(svm.SVC()))
decoder = SlidingEstimator(clf, scoring='roc_auc')
decoder.fit(X, y)
decoder
coef = get_coef(decoder, 'patterns_', inverse_transform=True)
coef = get_coef(decoder, 'filters_')
coef = get_coef(decoder, 'filters_', inverse_transform=True)
coef.shape

from mne import io, EvokedArray
# Extract and plot patterns and filters
for name in ('patterns_', 'filters_'):
    # The `inverse_transform` parameter will call this method on any estimator
    # contained in the pipeline, in reverse order.
    coef = get_coef(decoder, name, inverse_transform=True)
    evoked = EvokedArray(coef, epochs_subset.info, tmin=epochs_subset.tmin)
    evoked.plot_topomap(times=np.arange(-0.1, 0.8, 0.1), cmap='viridis')
    evoked.plot_topomap(times='interactive', cmap='viridis')




# evoked potentials (ERP) average
evoked_congruent = epochs_subset['congruent'].average()
evoked_incongruent = epochs_subset['incongruent'].average()

f, axs = plt.subplots(2, 1, figsize=(10, 5))
f.suptitle('congruent vs incongruent', fontsize=20)
evoked_congruent.plot(axes=axs[0], show=False, time_unit='s')
evoked_incongruent.plot(axes=axs[1], show=False, time_unit='s')
plt.tight_layout()
plt.show()

