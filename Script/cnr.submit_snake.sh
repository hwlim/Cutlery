#!/usr/bin/env bash

## Written by Hee Woong Lim and Chris Ahn
## Script to submit snakemake job to CCHMC lsf system

source $COMMON_LIB_BASE/commonBash.sh

nJob=50
totalWaitTime="48:00"
config=${CUTLERY}/Snakemake/cluster.yml

mkdir -p logs

## The job script is passed through stdin (here-document) instead of a bsub command string,
## since LSF expands $(...) in a command string before the job runs, i.e. before anaconda3 is loaded.
## Variables are filled in at submission; \$ is kept for the job to evaluate.
## <<- strips leading tabs, so the here-document must be indented with tabs only.

## if config.yml file exists in current directory, run default beginner mode
if [ -e "config.yml" ]; then
	echo -e "config.yml file found in work directory; running Cutlery in default mode." >&2

	bsub -W ${totalWaitTime} -q rhel9 -eo submit.err -oo submit.out <<- EOF
		module purge
		module load anaconda3/1.0.1
		eval "\$(conda shell.bash hook)"
		conda activate snakemake-7.18.2
		export XDG_CACHE_HOME=/scratch/$USER/snakemake-cache

		snakemake -j $nJob --rerun-incomplete \
			-s ${CUTLERY}/Snakemake/Snakefile_Default \
			--latency-wait 60 \
			--cluster-config $config \
			--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -q rhel9 -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'
	EOF

## run advanced mode if config.yml file doesn't exist in the current directory
else
	echo -e "No config.yml file found in current directory; running Cutlery in advanced mode." >&2

	bsub -W ${totalWaitTime} -q rhel9 -eo submit.err -oo submit.out <<- EOF
		module purge
		module load anaconda3/1.0.1
		eval "\$(conda shell.bash hook)"
		conda activate snakemake-7.18.2
		export XDG_CACHE_HOME=/scratch/$USER/snakemake-cache

		snakemake -j $nJob --rerun-incomplete \
			--latency-wait 60 \
			--cluster-config $config \
			--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -q rhel9 -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'
	EOF
fi
