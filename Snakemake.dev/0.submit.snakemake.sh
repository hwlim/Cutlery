#!/usr/bin/env bash

nJob=20
totalWaitTime="24:00"

snakemake --dag | dot -Tpdf > diag.pdf


mkdir -p logs
bsub -W ${totalWaitTime} -e submit.err -o submit.out \
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 30 \
		--cluster-config cluster.yml \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

