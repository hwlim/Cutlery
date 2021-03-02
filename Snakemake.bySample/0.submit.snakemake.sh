#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh


nJob=20
totalWaitTime="48:00"
timestamp=$(date +%Y%m%d_%H%M%S)
cluster_config=${CUTLERY}/Snakemake.bySample/cluster.yml

assertFileExist $cluster_config

#if [ ! -f diag.pdf ];then
if [  "`which python3`" == "" ]; then
	module load python3/3.6.3
else
	python_version=`python3 --version`
	[ "$python_version" != "Python 3.6.3" ] && module load python3/3.6.3
fi
module load graphviz/2.40.1
snakemake --dag | dot -Tpdf > diag.pdf
#fi

mkdir -p logs
bsub -W ${totalWaitTime} -eo bsub.err -oo bsub.out \
	"module load python3/3.6.3
	snakemake -j $nJob \
		--latency-wait 60 \
		--cluster-config $cluster_config \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'"

