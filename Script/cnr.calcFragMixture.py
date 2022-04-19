#!/usr/bin/env python3

# get frag QC info for sample and output as line separated by tabs
# nfrFrag%    nucFrag%    w0    w1      w2

import pandas as pd
import csv
import argparse
import gzip
import numpy as np
from scipy import stats
from scipy.stats import norm
from sklearn.mixture import GaussianMixture as GMM
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

options = argparse.ArgumentParser(description="Get fragment QC stats. This script uses a Gaussian mixture model to determine the weights of nfr fragments and nuc fragments.", usage="python3 getFragQC.py (options) fragLenDist.txt")
options.add_argument('fragLenDist_file',
                        help='Required; fragLenDist.txt file from Cutlery')
options.add_argument('-o', '--outputFile', default='fragMix',
                        help='Prefix for output files. There are three output files in total: output.txt file containing the ratios and weights of nfr and nuc fragments, along with a GMM.pdf and a GMM.png file visualizing the mixture model.')
args = options.parse_args()

# open density file as df
fragLenDist = pd.read_csv(args.fragLenDist_file, delimiter = "\t")

# get nfr and nuc fragment percentages
totalCount = fragLenDist['Cnt'].sum()
nfrCount = fragLenDist['Cnt'][:120].sum()
nucCount = fragLenDist['Cnt'][150:].sum()

# get fragment ratio
nfrFrac = nfrCount / totalCount * 100
nucFrac = nucCount / totalCount * 100

# open density file
fragLenDist = open(args.fragLenDist_file)
# convert density file to individual distances
densityFile = csv.reader(fragLenDist, delimiter="\t")
# skip header
next(densityFile)

tmpList = []
for row in densityFile:
    count = int(row[1])
    i = 0
    while i < count:
        tmpList.append(int(row[0]))
        i = i + 1

# convert to df
npArray = np.asarray(tmpList)
df = pd.DataFrame(npArray)
df = df[df[0]<1000].sample(1000000)

# Get weights of distribution from GMM
# define number of clusters/components as 3
comp = 3

# fit a Gaussian Mixture Model
gmm = GMM(n_components = comp, max_iter=100, random_state=10, covariance_type = 'full')
gmm.fit(df)

# Sort mean, covariances and weights
me = gmm.means_
mean = [me[0][0],me[1][0],me[2][0]]
mean = np.asarray(mean)
sortedInd = mean.argsort()

cov  = gmm.covariances_
covs = [cov[0][0][0],cov[1][0][0],cov[2][0][0]]
covs = np.asarray(covs)

weight = gmm.weights_
weights = [weight[0],weight[1],weight[2]]
weights = np.asarray(weights)

mean = mean[sortedInd[:]]
covs = covs[sortedInd[:]]
weights = weights[sortedInd[:]]

# Create and save fragQC table
outName = args.outputFile + ".txt"
with open(outName, mode = 'w') as outFile:
    outFile_write = csv.writer(outFile, delimiter='\t')
    outFile_write.writerow(["NFR_Frags(%)", "NUC_Frags(%)", "W_0", "W_1", "W_2"])
    outFile_write.writerow([nfrFrac, nucFrac, weights[0], weights[1], weights[2]])


#plot gmm
x_axis = np.arange(0, 1000, 2)
y_axis0 = norm.pdf(x_axis, float(mean[0]), np.sqrt(float(covs[0])))*weights[0] # 1st gaussian
y_axis1 = norm.pdf(x_axis, float(mean[1]), np.sqrt(float(covs[1])))*weights[1] # 2nd gaussian
y_axis2 = norm.pdf(x_axis, float(mean[2]), np.sqrt(float(covs[2])))*weights[2] # 3rd gaussian

plt.hist(df, density=True, color='black', bins=np.arange(0, 1000, 5))
plt.plot(x_axis, y_axis0, lw=3, c='C0', label = round(weights[0], 3))
plt.plot(x_axis, y_axis1, lw=3, c='C1', label = round(weights[1], 3))
plt.plot(x_axis, y_axis2, lw=3, c='C2', label = round(weights[2], 3))
#plt.plot(x_axis, y_axis0+y_axis1+y_axis2, lw=3, c='C3', ls='dashed', label = "Combined")
plt.legend()
#plt.xlim(0, 1000)
plt.xlabel(r"Distance", fontsize=20)
plt.ylabel(r"Density", fontsize=20)
plt.savefig(args.outputFile + ".GMM.png")
plt.savefig(args.outputFile + ".GMM.pdf")
