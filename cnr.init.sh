#!/usr/bin/env bash

if [ -z ${CUTLERY+x} ]; then
	echo -e "Error: Environment variable CUTLERY is not defined" >&2
	exit 1
fi

echo -e "Initializing Cutlery folder" >&2
echo -e " CUTLERY = ${CUTLERY}" >&2


cp -v ${CUTLERY}/Snakemake.bySample/sample.tsv .
cp -v ${CUTLERY}/Snakemake.bySample/Snakemake .
cp -v ${CUTLERY}/Snakemake.bySample/0.submit.snakemake.sh .
