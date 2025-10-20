#!/usr/bin/env bash
## Script to perform snakemake dry-run

if [ -z ${CUTLERY+x} ]; then
	echo -e "Error: Environment variable CUTLERY is not defined. Refer to the initial setup section in the github page at https://github.com/hwlim/Cutlery for instructions on setting up this enrironment variable." >&2
	exit 1
fi

echo -e "Performing dry-run..." >&2
module purge
module load python3/3.6.3

if [ -e "config.yml" ]; then
	snakemake -np -s ${CUTLERY}/Snakemake/Snakefile_Default
else
	snakemake -np
fi

module purge
