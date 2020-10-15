#!/usr/bin/env python
# coding: utf-8

# # Data exploration
# The notebook presents the preprocessed satellite data with general data exploration techniques.
# %matplotlib qt
import os
print(os.getcwd())
from os import path
import pandas
import numpy as np
from matplotlib import pyplot as plt
import random
import sys

sys.path.append(os.path.abspath("src/utils"))
os.environ["PROJ_LIB"] = "C:\\Utilities\\Python\\Anaconda\\Library\\share"; #fixr
import toolbox

try:
    os.mkdir('figures')
except:
    pass
# %% Loading preprocessed data
df = toolbox.load_data()


# %% Exploring fields
# Number of rows.
len(df)

# The preprocessing step already removed name fields from the NetCDF files. Only the fields useful for the study are present.
list(df.columns.values)

# Quick statistics about the data.
df.describe()

# The preprocessing removes all rows where the basics fields (time, lon/lat, tb, sig0_kua and snow_facies_flag) are set to NaN.
# But other fields might have invalid data.
df.isnull().sum(axis = 0)

# We are particularly interested in the snow_facies_flag.
# The cell below shows for which column data is available for each group.
df.groupby("ice_sheet_snow_facies_flag_01_ku").count()

# If waveform data is also available. Let's see some samples.
try:
    waveforms_ku = np.array(df[["waveform_20_ku_%d" % i for i in range(128)]].dropna())

    plt.figure(figsize=(10, 10))
    plt.title("Waveform samples")
    plt.plot(waveforms_ku[0])
    plt.plot(waveforms_ku[10000])
    plt.plot(waveforms_ku[20000])
    plt.plot(waveforms_ku[30000])
    plt.plot(waveforms_ku[40000])
    plt.show()
    
    # The mean wave
    plt.figure(figsize=(10, 10))
    plt.title("Mean wave")
    plt.plot(waveforms_ku.mean(axis=0))
    plt.show()
    
    # The samples show a varying amplitude. Let's plot the norm of the waveform as an histogram.
    plt.figure(figsize=(10, 10))
    plt.title("Histogram of waveform amplitude")
    plt.hist(np.linalg.norm(waveforms_ku, axis=1), range=(0, 12000),  bins=100)
    plt.show()
    
    # Waveforms have an amplitude that varies a lot. Plotting the normalized waveform will help gain insight into the shape of the waves.
    wave_norm = (waveforms_ku.T / np.linalg.norm(waveforms_ku, axis=1)).T
    wave_mean = wave_norm.mean(axis=0)
    wave_std = wave_norm.std(axis=0)
    plt.figure(figsize=(10, 10))
    plt.title("Mean and deviation of normalised waveforms.")
    plt.plot(wave_mean)
    plt.plot(wave_mean - wave_std, "b--")
    plt.plot(wave_mean + wave_std, "b--")
    plt.show()
    
    # Or as a boxplot:
    plt.figure(figsize=(10, 10))
    plt.title("Boxplot of normalized waveform")
    plt.boxplot(wave_norm, showfliers=False)
    plt.show()
except:
    pass

# %% Visualisations
#  Dataset TSNE
# T-stochastic Neighbourg Embedding is a useful tool to visualise high dimension data in 2D.
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

# Sampling a small sample for faster computing.
sample1000 = df.sample(n=1000)

# TSNE visualisation for the metrics from the original classification.
sample1000_filtered = sample1000[[ "tb_ratio", "tb_avg", "sig0_ice_sheet_01_ku", "ku_c_diff"]]
sample1000_filtered = (sample1000_filtered - sample1000_filtered.mean()) / sample1000_filtered.std()
data_emb = TSNE(random_state=42).fit_transform(sample1000_filtered)

fig, ax = plt.subplots(1,2, figsize=(12, 5))
ax[0].set_title("TSNE of Temperature ratio and average, \nsig ku and ku/c difference")
sc = ax[0].scatter(data_emb[:, 1], data_emb[:, 0], c=sample1000["ice_sheet_snow_facies_flag_01_ku"], cmap=plt.cm.plasma)
plt.colorbar(sc, ax=ax[0])

# TSNE visualisation for the original metrics and plrm waveform parameters. 
sample1000_filtered = sample1000[[
    "tb_ratio", "tb_avg", "sig0_ice_sheet_01_ku", "ku_c_diff", "sig0_leading_edge_ice_20_plrm_ku",
    "width_leading_edge_ice_20_plrm_ku", "slope_first_trailing_edge_ice_20_plrm_ku",
    "slope_second_trailing_edge_ice_20_plrm_ku", "peakiness_1_20_plrm_ku"
]].dropna()

sample1000_filtered = (sample1000_filtered - sample1000_filtered.mean()) / sample1000_filtered.std()
data_emb = TSNE(random_state=42).fit_transform(sample1000_filtered)

ax[1].set_title("TSNE of Temperature ratio and average, \nsig ku and ku/c difference and waveform parameters.")
sc = ax[1].scatter(data_emb[:, 1], data_emb[:, 0], c=sample1000["ice_sheet_snow_facies_flag_01_ku"].loc[sample1000_filtered.index.values], cmap=plt.cm.plasma)
plt.colorbar(sc, ax=ax[1])
fig.savefig('figures/TSNE.jpg', bbox_inches='tight',dpi=90)

# %%# Class caracterisation
# It is usefull to explore each classe for each parameter. The plots below show
# histograms for each field for each group and the full data.
grouped_by_class = df.groupby("ice_sheet_snow_facies_flag_01_ku")

def histbygroup(field, df_all, grouped, ax , range=None, bins=100):
    """Displays histogram for total data and data grouped per class for a given field."""
    ax.hist(df_all[field].dropna(), bins=bins, range=range)
    legends = ["All"]
    for name, group in grouped:
        ax.hist(group[field].dropna(), bins=bins, histtype="step", range=range)
        legends = ["group " + str(name)] + legends
    ax.legend(legends)
    ax.set_xlabel(field)

fig, ax =     plt.subplots(3,3, figsize=(15, 12))
fig.subplots_adjust(hspace = 0.22,wspace = 0.2)
histbygroup("tb_ratio", df, grouped_by_class, ax[0,0])
histbygroup("tb_avg",   df, grouped_by_class, ax[0,1])
histbygroup("sig0_ice_sheet_01_ku", df, grouped_by_class, ax[0,2])
histbygroup("ku_c_diff", df, grouped_by_class, ax[1,0])
histbygroup("sig0_leading_edge_ice_20_plrm_ku", df, grouped_by_class, ax[1,1])
histbygroup("width_leading_edge_ice_20_plrm_ku", df, grouped_by_class, ax[1,2], bins=30)
histbygroup("slope_first_trailing_edge_ice_20_plrm_ku", df, grouped_by_class, ax[2,0], range=(-5e7, 5e7))
histbygroup("slope_second_trailing_edge_ice_20_plrm_ku", df, grouped_by_class, ax[2,1], range=(-5e7, 5e7))
histbygroup("peakiness_1_20_plrm_ku", df, grouped_by_class, ax[2,2], range=(0, 20))
ax[1,0].set_ylabel('Counts')
#fig.tight_layout()
fig.savefig('figures/histograms_classes.jpg', bbox_inches='tight',dpi=120)


# %%# Geographic plots
# Geographic repartition of the classes.
def drawClasses(cycle):
    plt.title("Cycle %d: %s to %s" % (
        cycle,
        toolbox.ESA_time_to_datetime(df[df["cycle"] == cycle]["time_01"].min()).strftime("%d-%b-%Y"),
        toolbox.ESA_time_to_datetime(df[df["cycle"] == cycle]["time_01"].max()).strftime("%d-%b-%Y")
    ))
    toolbox.drawGreenland(df[df["cycle"] == cycle], "lon_01", "lat_01", "ice_sheet_snow_facies_flag_01_ku","Original classes")

fig = plt.figure(figsize=(12, 20))
plt.subplot(1, 2, 1)
drawClasses(40)
plt.subplot(1, 2, 2)
drawClasses(46)
plt.show()
fig.savefig('figures/map_classes.jpg', bbox_inches='tight',dpi=120)

#%% Geographic representation of the parameters.
fields_of_interest = [
    "tb_ratio", "tb_avg", "sig0_ice_sheet_01_ku", "ku_c_diff",
    "sig0_leading_edge_ice_20_plrm_ku", "width_leading_edge_ice_20_plrm_ku",
    "slope_first_trailing_edge_ice_20_plrm_ku", "slope_second_trailing_edge_ice_20_plrm_ku",
    "peakiness_1_20_plrm_ku"
]
simple_names = [
    "tb_ratio", "tb_avg", "sig0_ice_sheet_01_ku", "ku_c_diff",
    "sig0_leading_edge_ice_20_plrm_ku", "width_leading_edge_ice_20_plrm_ku",
    "slope_first_trailing_edge_ice_20_plrm_ku", "slope_second_trailing_edge_ice_20_plrm_ku",
    "peakiness_1_20_plrm_ku"
]

def drawParameters(cycle):
    fig = plt.figure(figsize=(13, 16))
    plt.suptitle("Cycle %d: %s to %s" % (
        cycle,
        toolbox.ESA_time_to_datetime(df[df["cycle"] == cycle]["time_01"].min()).strftime("%d-%b-%Y"),
        toolbox.ESA_time_to_datetime(df[df["cycle"] == cycle]["time_01"].max()).strftime("%d-%b-%Y")
    ))
    for i, field in enumerate(fields_of_interest):
        plt.subplot(3, 3, i + 1)
#        plt.title(field)
        toolbox.drawGreenland(df[df["cycle"] == cycle], "lon_01", "lat_01", field,simple_names[i])
    plt.show()
    fig.savefig('figures/map_param_cycle_'+str(cycle)+'.jpg', bbox_inches='tight',dpi=120)


drawParameters(40)
drawParameters(46)
# %% Correlation between parameters
# The map shows that *slope_first_trailing_edge_ice_20_plrm_ku* and 
# *slope_second_trailing_edge_ice_20_plrm_ku* seemed to be correlated, 
# as well as *sig0_leading_edge_ice_20_plrm_ku* and *sig0_ice_sheet_01_ku*, 
# which is unsuprising since they are actually very similar metrics.
# # Let's see how correlated they are.
print(df[["slope_first_trailing_edge_ice_20_plrm_ku", "slope_second_trailing_edge_ice_20_plrm_ku"]].corr())
# The trailing slopes are not actually correlated.
print(df[["sig0_leading_edge_ice_20_plrm_ku", "sig0_ice_sheet_01_ku"]].corr())
# But the sig0 are.
# Finally the full correlation matrix.
fields_of_interest = [
    "tb_ratio", "tb_avg", "sig0_ice_sheet_01_ku", "ku_c_diff", "sig0_leading_edge_ice_20_plrm_ku",
    "width_leading_edge_ice_20_plrm_ku", "slope_first_trailing_edge_ice_20_plrm_ku",
    "slope_second_trailing_edge_ice_20_plrm_ku", "peakiness_1_20_plrm_ku"
]
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
cax = ax.matshow(df[fields_of_interest].dropna().corr())
fig.colorbar(cax)

ax.set_xticklabels(['']+fields_of_interest, rotation=90)
ax.set_yticklabels(['']+fields_of_interest)
fig.savefig('figures/cor_matrix.jpg', bbox_inches='tight',dpi=120)



