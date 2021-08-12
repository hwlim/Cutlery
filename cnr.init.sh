#!/usr/bin/env bash

if [ -z ${COMMON_LIB_BASE+x} ]; then
	echo -e "Error: LimLabBase is not perperly set up" >&2
	exit 1
fi
source $COMMON_LIB_BASE/commonBash.sh

if [ -z ${CUTLERY+x} ]; then
	echo -e "Error: Environment variable CUTLERY is not defined" >&2
	exit 1
fi

echo -e "Initializing Cutlery analysis" >&2
echo -e "  CUTLERY path = ${CUTLERY}" >&2


cp -v ${CUTLERY}/Snakemake.bySample/sample.tsv .
cp -v ${CUTLERY}/Snakemake.bySample/Snakefile .
#cp -v ${CUTLERY}/Snakemake.bySample/0.submit.snakemake.sh .
