#!/usr/bin/env bash

## Written by Hee Woong Lim
##
## Script to submit snakemake job to CCHMC lsf system

source $COMMON_LIB_BASE/commonBash.sh


nJob=30
totalWaitTime="48:00"
#timestamp=$(date +%Y%m%d_%H%M%S)
config=${CUTLERY}/Snakemake/cluster.yml

assertFileExist $config
assertFileExist ./Snakefile



if [ ! -f diag.pdf ];then
	module load python3/3.6.3
	module load graphviz/2.40.1
	snakemake --dag | dot -Tpdf > diag.pdf
fi

#module load python3/3.6.3
#snakemake -np
#exit 0
mkdir -p logs
bsub -W ${totalWaitTime} -eo submit.err -oo submit.out \
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 60 \
		--cluster-config $config \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

