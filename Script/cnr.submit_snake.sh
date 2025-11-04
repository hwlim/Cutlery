#!/usr/bin/env bash

## Written by Hee Woong Lim and Chris Ahn
## Script to submit snakemake job to CCHMC lsf system

source $COMMON_LIB_BASE/commonBash.sh

nJob=50
totalWaitTime="48:00"
config=${CUTLERY}/Snakemake/cluster.yml

mkdir -p logs

## if config.yml file exists in current directory, run default beginner mode
if [ -e "config.yml" ]; then
	echo -e "config.yml file found in work directory; running Cutlery in default mode." >&2
	
	if [ ! -f diag.pdf ];then
	module load python3/3.6.3
	module load graphviz/2.40.1
	snakemake -s ${CUTLERY}/Snakemake/Snakefile_Default --dag | dot -Tpdf > diag.pdf
	module purge
	fi
	
	bsub -W ${totalWaitTime} -q rhel9 -eo submit.err -oo submit.out \
		"
		module purge
		module load anaconda3
		source activate snakemake-7.18.2
		module load squashfs-tools/4.5.0
		snakemake -j $nJob --rerun-incomplete \
			-s ${CUTLERY}/Snakemake/Snakefile_Default \
			--latency-wait 60 \
			--cluster-config $config \
			--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -q rhel9 -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'
		"

## run advanced mode if config.yml file doesn't exist in the current directory
else
	echo -e "No config.yml file found in current directory; running Cutlery in advanced mode." >&2
	
	if [ ! -f diag.pdf ];then
	module load python3
	module load graphviz/2.40.1
	snakemake --dag | dot -Tpdf > diag.pdf
	fi
	
	bsub -W ${totalWaitTime} -q rhel9 -eo submit.err -oo submit.out \
		"
		module purge
		module load anaconda3
		source activate snakemake-7.18.2
		module load squashfs-tools/4.5.0
		snakemake -j $nJob \
			--latency-wait 60 \
			--cluster-config $config \
			--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -q rhel9 -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'
		"
fi
