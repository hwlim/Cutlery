#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh


nJob=20
totalWaitTime="48:00"
timestamp=$(date +%Y%m%d_%H%M%S)
cluster_config=~/bin/CnR/Snakemake.bySample/cluster.yml
src_snakefile=~/bin/CnR/Snakemake.bySample/Snakefile

if [ ! -f diag.pdf ];then
	module load python3/3.6.3
	snakemake --dag --snakefile $src_snakefile | dot -Tpdf > diag.pdf
fi

mkdir -p logs
bsub -W ${totalWaitTime} -eo bsub.${timestamp}.err -oo bsub.${timestamp}.out \
	"module load python3/3.6.3
	snakemake --snakefile $src_snakefile \
		-j $nJob \
		--latency-wait 60 \
		--cluster-config $cluster_config \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

