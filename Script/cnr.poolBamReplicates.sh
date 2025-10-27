#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Pool bam files of multiple replicates by group
#
# NOTE:
# Originally developed for CUT&Run pipeline. Cutlery
# But, now Cutlery uses fragment files, thus has been integrated to ngs.poolBamReplicates.sh with -s option
# Will be deprecated


source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <sample.tsv> <src bam directory> <des bam directory>
Description:
	Merge multiple bam files of a group according to sample/group information within a given sample.tsv file
Input:
	- sample.tsv file: containing columns 'Name' and 'Group', from original/replicate Cutlery run
	- src bam directory: contaning replicate bam files
	  e.g. <src sample dir>
	  ├── <sample1>.bam
	  └── <sample2>.bam
	- des bam directory: to write merged bam files
	  e.g. <des group dir>
	  ├── <group1>.bam
	  └── <group2>.bam
Options:
	-t: if set, dry run simply displaying pooling message, default=off
	-b: if set, bsub are submitted for merging bam files, default=off
	-u: if set, assuming unsorted replicate bam file (name-sorted as in old version of Cutlery), default=off
		FYI, previosu Cutlery created unsorted bam file.
		But the current version create coordinate-sorted bam files after alignment
	-f: if set, force overwrite existing destination bam files, default=off" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
testOnly=FALSE
bsub=FALSE
unsortedBam=FALSE
overwrite=FALSE
while getopts ":buft" opt; do
	case $opt in
		t)
			testOnly=TRUE
			;;
		b)
			bsub=TRUE
			;;
		u)
			unsortedBam=TRUE
			;;
		f)
			overwrite=TRUE
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
if [ $# -lt 3 ];then
	printUsage
	exit 1
fi
sampleInfo=$1
srcDir=$2
desDir=$3

assertFileExist $sampleInfo
assertDirExist $srcDir

###################################
## main code

#Usage: ngs.concateBamFiles.sh (options) [bam1] ...
#Description:
#	Concatenate bam files. No sorting
#Options:
#	-o <out>: Output file prefix including pathr. required
#	-m <mem>: Memory for sorting. default=5G
#	-s : If set, coordinate-sorted and indexed

groupL=`tail -n +2 $sampleInfo | grep -v -e ^$ -e ^# | cut -f 3 | sort | uniq`

echo -e "Pooling replicate bam files" >&2
echo -e "  - sampleInfo: $sampleInfo" >&2
echo -e "  - srcDir: $srcDir" >&2
echo -e "  - desDir: $desDir" >&2
echo -e "  - testOnly: $testOnly" >&2
echo -e "  - unsorted bam: $unsortedBam" >&2
echo -e "  - bsub:   $bsub" >&2


for group in ${groupL[@]}
do
	mkdir -p $desDir/${group}
	des=${desDir}/${group}/align.bam
	log=${desDir}/${group}/align.log

	## Checking existing destination file
	if [ -f $des ] && [ "$overwrite" != "TRUE" ];then
		echo -e "Warning: $des already exists, pass" >&2
		continue
	fi

	## List of replicate bam files
	srcL=( `tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s/align.bam\n", $2 }'` )
	assertFileExist ${srcL[@]}

	## Merging to destination
	echo -e "Creating $des" >&2
	for src in ${srcL[@]}
	do
		echo -e "  - $src" >&2
	done 2>&1 | tee $log

	if [ "$testOnly" == "TRUE" ];then
		continue
	fi

	if [ "$unsortedBam" == "TRUE" ];then
		optStr=""
	else
		optStr="-S"
	fi

	if [ "$bsub" == "TRUE" ];then
		## Parallel processing using HPC:lsf
		bsub -W 24:00 -n 1 "module load samtools/1.9.0; ngs.concateBamFiles.sh $optStr -o $des ${srcL[@]}"
	else
		## Sequential processing
		ngs.concateBamFiles.sh $optStr -o $des ${srcL[@]}
	fi
done

