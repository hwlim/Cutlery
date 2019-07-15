#!/usr/bin/env bash

#########################$$$$$$$$$$$$$$$
## Written by Hee-Woong Lim
#
## Wrapper script for star alignment
##
## To do:
## - output prefix or output file? for convenience in snakemake

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/_temp_.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/_temp_.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [fq1] (fq2)" >&2
	echo "Description: Wrapper script for star alignnment
Options:
	-g <reference>: STAR index directory (required)
	-o: outPrefix including path, default=align
		Since STAR does not add dot (.) after prefix, this script automatically add one if it doesn't end with dot
	-t: number of threads for parallel processing, default=4
	-s: sort by coordinate, defalt=Off
	-p: option string for STAR, default=NULL
Option examples:
	ChIP-seq / ATAC-seq:
		--alignSJDBoverhangMin 999 --alignIntronMax 1 --alignMatesGapMax 2000
		--outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.05 --outReadsUnmapped None
	Additional option to prevent soft-clipping:
		--alignEndsType EndToEnd 
	RNA-seq:
		--outSAMstrandField intronMotif --outFilterMultimapNmax 1 --outReadsUnmapped None
" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handle
genome=NULL
outPrefix=STAR
thread=4
sortBam=FALSE
optStr=NULL
while getopts ":g:o:t:p:s" opt; do
	case $opt in
		g)
			genome=$OPTARG
			;;
		o)
			outPrefix=$OPTARG
			;;
		t)
			thread=$OPTARG
			;;
		s)
			sortBam=TRUE
			;;
		p)
			optStr=$OPTARG
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

###########################
### Option handling
shift $((OPTIND-1))
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

if [ $# -eq 1 ];then
	fq1=$1
	fq2=NULL
	assertFileExist $fq1
else
	fq1=$1
	fq2=$2
	assertFileExist $fq1 $fq2
fi

if [ "$genome" = "NULL" ];then
	echo "Error: STAR reference directory must be specified (-g)" >&2
	printUsage
	exit 1
fi
isDirExist $genome

if [ "$sortBam" == "TRUE" ];then
	optStr="${optStr} --outSAMtype BAM SortedByCoordinate"
else
	optStr="${optStr} --outSAMtype BAM Unsorted"
fi


[[ ! "$outPrefix" =~ \.$ ]] && outPrefix=${outPrefix}.
desDir=`dirname $outPrefix`
mkdir -p $desDir
desBam=${outPrefix}bam


###########################
### Start processing

echo -e "============================================"
echo -e "Running STAR" >&2
echo -e "- fq1: $fq1" >&2
echo -e "- fq2: $fq2" >&2
echo -e "- Index: $genome" >&2
echo -e "- sortBam: $sortBam" >&2
echo -e "- outPrefix: $outPrefix" >&2
echo -e "- desBam: $desBam" >&2


if [ -f ${desBam} ];then
	echo -e "Error: ${desBam} already exists" >&2
	exit 1
fi


if [ "${fq2}" == "NULL" ];then
	echo -e "STAR --runMode alignReads
	--genomeDir ${genome}
	--readFilesIn <( zcat $fq1 )
	--genomeLoad NoSharedMemory
	--outFileNamePrefix ${outPrefix}
	${optStr}" >&2
	STAR --runMode alignReads \
		--genomeDir ${genome} \
		--readFilesIn <( zcat $fq1 ) \
		--genomeLoad NoSharedMemory \
		--outFileNamePrefix ${outPrefix} \
		${optStr}
else
	echo -e "STAR --runMode alignReads
	--genomeDir ${genome}
	--readFilesIn <( zcat $fq1 ) <( zcat $fq2 )
	--genomeLoad NoSharedMemory
	--outFileNamePrefix ${outPrefix}
	${optStr}" >&2
	STAR --runMode alignReads \
		--genomeDir ${genome} \
		--readFilesIn <( zcat $fq1 ) <( zcat $fq2 ) \
		--genomeLoad NoSharedMemory \
		--outFileNamePrefix ${outPrefix} \
		${optStr}
fi

if [ "$sortBam" = "TRUE" ];then
	mv ${outPrefix}Aligned.sortedByCoord.out.bam ${desBam}
	echo -e "Indexing" >&2
	samtools index ${desBam}
else
	mv ${outPrefix}Aligned.out.bam ${desBam}
fi

