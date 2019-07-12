#!/usr/bin/env bash

nthread=8

bsub -W 6:00 -n $nthread -e submit.err -o submit.out \
	"module load python3/3.6.3
	snakemake -j $nthread --latency-wait 30 --cluster-config cluster.json  --cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J {cluster.name} -eo {cluster.error} -oo {cluster.output}'"

