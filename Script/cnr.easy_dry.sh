#!/usr/bin/env bash

## Written by Hee Woong Lim
##
## Script to submit snakemake job to CCHMC lsf system
## Easy version 

source $COMMON_LIB_BASE/commonBash.sh

cutlery_config=./config.yml
cluster_config=${CUTLERY}/Snakemake/cluster.yml
snakefile=${CUTLERY}/Snakemake/Snakefile_Easy

assertFileExist $cutlery_config
assertFileExist $cluster_config
assertFileExist $snakefile

module purge
module load python3/3.6.3
snakemake -np --snakefile $snakefile

