#!/usr/bin/env bash

## Written by Hee Woong Lim
##
## Script to submit snakemake job to CCHMC lsf system
## Easy version 

source $COMMON_LIB_BASE/commonBash.sh

nJob=50
totalWaitTime="48:00"

cutlery_config=./config.yml
cluster_config=${CUTLERY}/Snakemake/cluster.yml
snakefile=${CUTLERY}/Snakemake/Snakefile_Easy

assertFileExist $cutlery_config
assertFileExist $cluster_config
assertFileExist $snakefile

mkdir -p logs
bsub -W ${totalWaitTime} -q normal -M 1000 -eo cutlery.err -oo cutlery.out <<- EOF
		module purge
		module load python3/3.6.3
		snakemake -j $nJob \
				--snakefile $snakefile \
				--latency-wait 60 \
				--cluster-config $cluster_config \
				--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -q normal -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'

EOF
