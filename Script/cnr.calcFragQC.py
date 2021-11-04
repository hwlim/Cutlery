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
# import astropy
# import matplotlib as mpl

options = argparse.ArgumentParser(description="Get fragment QC stats", usage="python3 getFragQC.py (options)")
options.add_argument('-f', '--fragmentBed_file',
                        help='frag.all.bed file')
options.add_argument('-d', '--fragLenDist_file',
                        help='fragLenDist.txt file from Cutlery')
options.add_argument('-o', '--outputFile',
                        help='Path to fragQC file')
args = options.parse_args()


# get nfr and nuc fragment percentages
totalCount = 0
nfrCount = 0
nucCount = 0

with gzip.open(args.fragmentBed_file, mode="rt") as fragBed:
    for frag in fragBed:
        
        totalCount = totalCount + 1

        splitStr = frag.split("\t")
        start = int(splitStr[1])
        end = int(splitStr[2])
        fragLen = end - start
        
        if fragLen <= 120:
            nfrCount = nfrCount + 1
        elif fragLen >= 150:
            nucCount = nucCount + 1

nfrFrac = nfrCount / totalCount * 100
nucFrac = nucCount / totalCount * 100

# convert density file to individual distances
fragLenDist = open(args.fragLenDist_file)
densityFile = csv.reader(fragLenDist, delimiter="\t")
# skip header
next(densityFile)

tmpList = []

check = 0

for row in densityFile:
    count = int(row[1])
    i = 0
    check = check + 1
    while i < count:
        tmpList.append(int(row[0]))
        #print(row[0], temp)
        i = i + 1

npArray = np.asarray(tmpList)

df = pd.DataFrame(npArray)

# Get weights of distribution from GMM
# define number of clusters/components
comp = 3

# fit a Gaussian Mixture Model with two components
gmm = GMM(n_components = comp, max_iter=100, random_state=10, covariance_type = 'full')
gmm.fit(df)

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

outName = args.outputFile + "/frag.QC.txt"
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
plt.savefig(args.outputFile + "/GMM.png")
plt.savefig(args.outputFile + "/GMM.pdf")
