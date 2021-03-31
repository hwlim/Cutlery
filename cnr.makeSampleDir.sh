#!/usr/bin/env bash

## Convert existing Cutlery processing results into the current version of sample folders
## Organizing all the output files in sample by sample manner not by output type 
## To save storage, existing out files/folders do not change. Instead, symbolic links are generated

source $COMMON_LIB_BASE/commonBash.sh


getDir()
{
	# $1 : Snakefile
	# $2 : variable name to retrieve
	src=$1
	name=$2

	set +o pipefail
	result=`grep ${name} $src | head -n 1 | sed -e 's/ //g' -e 's/=/\t/' | cut -f 2`
	set -o pipefail

	echo $result
}


nameL=(
bigWigDir1bp
bigWigDir1bp_abs
bigWigDirAllFrag
bigWigDir
filteredDir
homerDir
)

if [ $# -lt 2 ];then
	echo -e "Usage: `basename $0` <sample.tsv> <Snakefile>" >&2
	exit 1
fi

src_sample=$1
src_smk=$2
assertFileExist $src_sample
assertFileExist $src_smk

sampleL=( `tail -n +2 $src_sample | grep -v -e "^#" -e "^$" | cut -f 2` )

for sample in ${sampleL[@]}
do
	echo $sample
done



