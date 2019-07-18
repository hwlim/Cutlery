#!/usr/bin/env bash

nJob=8
totalWallTime="24:00"


mkdir -p logs
#bsub -W 6:00 -n $nthread -e submit.err -o submit.out \
bsub -W ${totalWallTime} -eo submit.err -oo submit.out \
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 30 \
		--cluster-config cluster.yml \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

