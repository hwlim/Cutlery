#!/usr/bin/env bash

nJob=8


#bsub -W 6:00 -n $nthread -e submit.err -o submit.out \
bsub -W 6:00 -e submit.err -o submit.out \ # no need to setup this because individual job will have their own threads
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 30 \
		--cluster-config cluster.yml \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J {cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

