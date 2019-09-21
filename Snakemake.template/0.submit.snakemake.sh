#!/usr/bin/env bash

nJob=20
totalWaitTime="48:00"

if [ ! -f diag.pdf ];then
	module load python3/3.6.3
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
		--cluster-config cluster.yml \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

