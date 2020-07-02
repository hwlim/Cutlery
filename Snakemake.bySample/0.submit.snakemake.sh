#!/usr/bin/env bash

nJob=20
totalWaitTime="48:00"
timestamp=$(date +%Y%m%d_%H%M%S)
cluster_config=~/bin/CnR/Snakemake.bySample/cluster.yml

if [ ! -f diag.pdf ];then
	module load python3/3.6.3
	snakemake --dag | dot -Tpdf > diag.pdf
fi

#module load python3/3.6.3
#snakemake -np
#exit 0
mkdir -p logs
bsub -W ${totalWaitTime} -eo bsub.${timestamp}.err -oo bsub.${timestamp}.out \
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 60 \
		--cluster-config $cluster_config \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

