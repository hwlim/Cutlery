#!/usr/bin/env bash
## Script to perform snakemake dry-run

if [ -z ${CUTLERY+x} ]; then
	echo -e "Error: Environment variable CUTLERY is not defined. Refer to the initial setup section in the github page at https://github.com/hwlim/Cutlery for instructions on setting up this enrironment variable." >&2
	exit 1
fi

module purge
module load anaconda3
source activate snakemake-7.18.2
module load squashfs-tools/4.5.0

if [ -e "config.yml" ]; then
	echo -e "Performing dry-run in Cutlery default mode..." >&2
	snakemake -np -s ${CUTLERY}/Snakemake/Snakefile_Default
else
	echo -e "Performing dry-run in Cutlery advanced mode..." >&2
	snakemake -np
fi

module purge
