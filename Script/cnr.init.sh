#!/usr/bin/env bash

## check environment path
if [ -z ${COMMON_LIB_BASE+x} ]; then
	echo -e "Error: LimLabBase is not properly set up. Check https://github.com/hwlim/Cutlery for setup information." >&2
	exit 1
fi
source $COMMON_LIB_BASE/commonBash.sh

## check for flag; -a flag is advanced mode and will create a snakemake file in the directory
if [[ "${1:-}" == "-a" ]]; then
	echo -e "Initializing advanced Cutlery." >&2
	echo -e "  CUTLERY path = ${CUTLERY}" >&2
	cp -i -v ${CUTLERY}/Snakemake/sample.tsv .
	cp -i -v ${CUTLERY}/Snakemake/Snakefile .

## default beginner mode will create config.cnr.yml file and no snakemake file
else
    echo -e "Initializing default Cutlery. Use flag -a for advanced Cutlery." >&2
	echo -e "  CUTLERY path = ${CUTLERY}" >&2
    cp -i -v ${CUTLERY}/Snakemake/sample.tsv .
	cp -i -v ${CUTLERY}/Snakemake/config.yml .
fi
