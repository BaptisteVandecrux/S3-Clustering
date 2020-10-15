#!/usr/bin/env python
# coding: utf-8

# # Clustering
# This notebook implements the clustering method described in "Snow Facies Over 
# Ice Sheets Derived From Envisat Active and Passive Observations
# https://www.doi.org/10.1109/TGRS.2008.2000818 and improves upon it.

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys
import skfuzzy as fuzz
from sklearn.metrics import confusion_matrix
sys.path.append(os.path.abspath("src/utils"))
import toolbox

df = toolbox.load_data()
# %% Original Clustering
# The original article gives the cluster centroids, but also the mean and 
# standard deviation used. The classificaiton can be done according to them. 

# Centroids for the 6 original groups
# x1: sig0_ku, x2: tb_avg, x3: tb_ratio, x4: sig0 diff ku/S 
envisat_cntr = np.array([
    [-1.7687, 1.6004, 1.0223, -1.6125],
    [-0.0572, -0.9718, -1.1724, 0.0599],
    [-1.0752, 0.1278, 0.2710, -0.8156],
    [0.3435, 0.5521, 0.4000, 0.1178],
    [0.8021, -0.3529, 0.3066, 0.6574],
    [0.2705, -0.0926, -0.3545, 0.1738]
])

envisat_means = np.array([ [8.5433, 190.0169, -0.0147, -0.6846] ])

envisat_stds = np.array([ [3.8460, 23.0642, 0.0153, 2.3071] ])

# Since the original article uses the S band of Envisat but Sentinel-3 has a 
# C-band instead, 2 clustering are done. One by ignoring the *sig0* different, 
# one using the Ku/C difference instead.

# Fuzzy C-means parameters
m = 2
err = 0.005
maxiter = 1000
ncenters = 6

# Preparing numpy data
data = np.vstack((df["sig0_ice_sheet_01_ku"], df["tb_avg"], df["tb_ratio"], df["ku_c_diff"]))

# Normalising
data = ((data.T - envisat_means) / envisat_stds).T

# With all 4 variables
u, u0, d, jm, p, fpc = fuzz.cluster.cmeans_predict(data, envisat_cntr, m, 
                                                   error=err, maxiter=maxiter, 
                                                   init=None, seed=42)
cluster_membership_4vars = np.argmax(u, axis=0)

# With 3 variables
u, u0, d, jm, p, fpc = fuzz.cluster.cmeans_predict(data[:3, :], envisat_cntr[:, :3], m,
                                                   error=err, maxiter=maxiter, 
                                                   init=None, seed=42)
cluster_membership_3vars = np.argmax(u, axis=0)

#Let's compare those classification to the provided classification.

fig = plt.figure(figsize=(8,8))
mat = confusion_matrix(cluster_membership_4vars, df["ice_sheet_snow_facies_flag_01_ku"] - 1, 
                       normalize="true")
plt.title("Confusion matrix between S3 flags and new classification (4 vars,  normalized from Envisat mean and std)")
plt.xlabel('ice_sheet_snow_facies_flag_01_ku - 1')
plt.ylabel('new classes (4 vars)')
plt.imshow(mat)
fig.savefig('figures/cor_matrix_classification_4var.jpg', bbox_inches='tight',dpi=120)

mat = confusion_matrix(cluster_membership_3vars, df["ice_sheet_snow_facies_flag_01_ku"] - 1, normalize="true")
fig = plt.figure(figsize=(8,8))
plt.title("Confusion matrix between S3 flags and new classification (3 vars,  normalized from Envisat mean and std)")
plt.xlabel('ice_sheet_snow_facies_flag_01_ku - 1')
plt.ylabel('new classes (3 vars)')
plt.imshow(mat)
fig.savefig('figures/cor_matrix_classification_3var.jpg', bbox_inches='tight',dpi=120)


# The classification does not seem to work properly. This can be due to a 
# difference in equipement and hence in the mean value and deviation of the measurements.

# %% This problem can be fixed by using computed mean and standard deviations 
# instead of the one provided for Envisat. 

data = np.vstack((df["sig0_ice_sheet_01_ku"], df["tb_avg"], df["tb_ratio"], df["ku_c_diff"]))

# Normalising
means = data.mean(axis=1)
stds = data.std(axis=1)
data = ((data.T - means) / stds).T

# With all 4 variables
u, u0, d, jm, p, fpc = fuzz.cluster.cmeans_predict(data, envisat_cntr, m, error=err, maxiter=maxiter, init=None, seed=42)
cluster_membership_4vars = np.argmax(u, axis=0)

# With 3 variables
u, u0, d, jm, p, fpc = fuzz.cluster.cmeans_predict(data[:3, :], envisat_cntr[:, :3], m, error=err, maxiter=maxiter, init=None, seed=42)
cluster_membership_3vars = np.argmax(u, axis=0)

mat = confusion_matrix(cluster_membership_4vars, df["ice_sheet_snow_facies_flag_01_ku"] - 1, normalize="true")
fig = plt.figure(figsize=(8,8))
plt.title("Confusion matrix between S3 flags and new classification (4 vars, normalized)")
plt.xlabel('ice_sheet_snow_facies_flag_01_ku - 1')
plt.ylabel('new classes (4 vars)')
plt.imshow(mat)
fig.savefig('figures/cor_matrix_classification_4var_norm.jpg', bbox_inches='tight',dpi=120)

mat = confusion_matrix(cluster_membership_3vars, df["ice_sheet_snow_facies_flag_01_ku"] - 1, normalize="true")
fig = plt.figure(figsize=(8,8))
plt.title("Confusion matrix between S3 flags and new classification (3 vars, nomalized)")
plt.xlabel('ice_sheet_snow_facies_flag_01_ku - 1')
plt.ylabel('new classes (3 vars)')
plt.imshow(mat)
fig.savefig('figures/cor_matrix_classification_3var_norm.jpg', bbox_inches='tight',dpi=120)


# Now the classification is very similar to the flag provided by ESA.
# The small difference can be explained by the difference in normalisation in 
# they used a differente timeframe for the mean and std. Unfortunatelly, the 
# technical documentation does not seem to provide the mean and standard 
# deviation used to compute the flag (nor does it confirms that the centroids 
# used are indeed the one provided in the paper).

# %% Reimplementation
# Since the classification is using different instruments than the ones on 
# Envisat, it might be useful to recompute the centroids for the Sentinel data.
# We'll be using the year 2017 as training data.

from importlib import reload
toolbox = reload(toolbox)

df2017 = df[
    (df["time_01"] >= toolbox.str_to_ESA_time("2017-01-01 00:00:00"))
    & (df["time_01"] < toolbox.str_to_ESA_time("2018-01-01 00:00:00"))
]

df2017.count()

# Preparing numpy data
data = np.vstack((df2017["sig0_ice_sheet_01_ku"], df2017["tb_avg"], df2017["tb_ratio"], df2017["ku_c_diff"]))

# Normalising
means2017 = data.mean(axis=1)
stds2017 = data.std(axis=1)
data = ((data.T - means2017) / stds2017).T

cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(data, ncenters, m, error=err, maxiter=maxiter, init=None, seed=42)
cluster_membership_reimpl = np.argmax(u, axis=0)

cntr

mat = confusion_matrix(cluster_membership_reimpl, df2017["ice_sheet_snow_facies_flag_01_ku"] - 1, normalize="true")
fig = plt.figure(figsize=(8,8))
plt.title("Confusion matrix between S3 flags and new classification (reimplementation)")
plt.xlabel('ice_sheet_snow_facies_flag_01_ku - 1')
plt.ylabel('new classes (reimplementation)')
plt.imshow(mat)
fig.savefig('figures/cor_matrix_classification_reimplemented.jpg', bbox_inches='tight',dpi=120)


df2017.groupby("cycle")["time_01"].count()


# %% mapping reimplementation
def drawCycle(dataframe, cycle, cluster_membership_reimpl,tag):
    fig = plt.figure(figsize=(10, 8))
    plt.suptitle("Map of classes on cycle %d. (from %s to %s)" % (
        cycle,
        toolbox.ESA_time_to_datetime(dataframe[dataframe["cycle"] == cycle]["time_01"].min()).strftime("%b %Y"),
        toolbox.ESA_time_to_datetime(dataframe[dataframe["cycle"] == cycle]["time_01"].max()).strftime("%b %Y")
    ))
    plt.subplot(1, 2, 1)
    plt.title("Reimplementation")
    toolbox.drawGreenland(dataframe[dataframe["cycle"] == cycle], "lon_01", "lat_01", 
                          cluster_membership_reimpl[dataframe["cycle"] == cycle],
                          'Reimplemented classes')
    plt.subplot(1, 2, 2)
    plt.title("Old")
    toolbox.drawGreenland(dataframe[dataframe["cycle"] == cycle], "lon_01", 
                          "lat_01", "ice_sheet_snow_facies_flag_01_ku",
                          'Original classes')
    fig.savefig('figures/map_classification_'+tag+'_'+str(cycle)+'.jpg', bbox_inches='tight',dpi=120)

drawCycle(df2017, 13, cluster_membership_reimpl,'reimplemented')
drawCycle(df2017, 20, cluster_membership_reimpl,'reimplemented')


# %% Improvements: Include waveform data
# The original classification does not contain waveform caracteristics. 
# To improve upon the classification, waveform parameters are added to the 
# fuzzy clustering.
# 
# Namely, *sig0* and *width* of the leading edge, slopes of the trailing edge 
# and peakiness of the wave.

# Preparing the data, filtering NaN
df2017_filtered = df2017[[
    "time_01",
    "lon_01",
    "lat_01",
    "cycle",
    "sig0_ice_sheet_01_ku",
    "tb_avg",
    "tb_ratio",
    "ku_c_diff",
    "sig0_leading_edge_ice_20_plrm_ku",
    "width_leading_edge_ice_20_plrm_ku",
    "slope_first_trailing_edge_ice_20_plrm_ku",
    "slope_second_trailing_edge_ice_20_plrm_ku",
    "peakiness_1_20_plrm_ku",
    "ice_sheet_snow_facies_flag_01_ku"
]].dropna()
df2017_filtered.count()

# Preparing numpy data
data = np.array(df2017_filtered[[
    "sig0_ice_sheet_01_ku",
    "tb_avg",
    "tb_ratio",
    "ku_c_diff",
    "sig0_leading_edge_ice_20_plrm_ku",
    "width_leading_edge_ice_20_plrm_ku",
    "slope_first_trailing_edge_ice_20_plrm_ku",
    "slope_second_trailing_edge_ice_20_plrm_ku",
    "peakiness_1_20_plrm_ku"
]]).T

# Normalising
means2017 = data.mean(axis=1)
stds2017 = data.std(axis=1)
data = ((data.T - means2017) / stds2017).T

# Clustering
cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(data, ncenters, m, error=err, maxiter=maxiter, init=None, seed=42)
cluster_membership_improved = np.argmax(u, axis=0)

means2017
stds2017
cntr

# %% Confusion matrix
mat = confusion_matrix(cluster_membership_improved, df2017_filtered["ice_sheet_snow_facies_flag_01_ku"] - 1, normalize="true")
plt.title("Confusion matrix between S3 flags and new classification (with waveform data)")
plt.xlabel('ice_sheet_snow_facies_flag_01_ku - 1')
plt.ylabel('new classes (with waveform data)')
plt.imshow(mat)
plt.show()


# %% maps
drawCycle(df2017_filtered, 13, cluster_membership_improved,'improved')
drawCycle(df2017_filtered, 20, cluster_membership_improved,'improved')
