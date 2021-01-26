#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Pool bam files of multiple replicates by group

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <sample.tsv> <src bam directory> <des bam directory>
Description:
	Merge multiple bam files of a group according to sample/group information within a given sample.tsv file
Input:
	- sample.tsv file: containing columns 'Name' and 'Group'
	- src bam directory: contanin replicate bam files
	- des bam directory: to write merged bam files
		bam files are named as <group>.bam
Options:
	-b: if set, bam files are merged in parallel by submitting multiple bsub jobs, default=off" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
bsub=0
while getopts ":b" opt; do
	case $opt in
		b)
			bsub=1
			;;
		\?)
			echo "Invalid options: -$OPTARG" >&2
			printUsage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			printUsage
			exit 1
			;;
	esac
done


shift $((OPTIND-1))
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi
sampleInfo=$1
srcDir=$2
desDir=$3


###################################
## main code

Usage: ngs.concateBamFiles.sh (options) [bam1] ...
Description:
	Concatenate bam files. No sorting
Options:
	-o <out>: Output file prefix including pathr. required
	-m <mem>: Memory for sorting. default=5G
	-s : If set, coordinate-sorted and indexed


