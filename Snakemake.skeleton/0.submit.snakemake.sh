#!/usr/bin/env bash

#module load python3/3.6.3

#snakemake -j 100 --cluster-config cluster.json --cluster "bsub -W {cluster.time} M {cluster.mem} -n {threads} -R span[ptile=4]"
bsub -W 1:00 "module load python3/3.6.3; snakemake -j 100 --cluster 'bsub' --latency-wait 30"

