#!/usr/bin/env bash

# ATAC-seq tools
# Written by Hee-Wooong Lim
# 
# Deduplication of ATAC-seq bam file using Picard tools
#	Sort by coordinate -> deduplicatation -> Sort by read names
#
# - Input: BAM, not necessarily sorted
# - Output: Produce two bam files
#	    1) Coordiated sorted bam file *.csort.bam and index
#	    2) Deduplicated bam file, *dedup.bam sorted by read name
#	    Duplicate read-pairs are deleted explicitly or simply marked as 'duplicate' in sam flag depending on the option
#



source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [inputFiles] ..." >&2
	echo -e "Description: Deduplicate given bam files using Picard tool" >&2
	echo -e "\t1) Coordinate sorting => 2) Deduplication => 3) Readname sorting" >&2
	echo -e "\tDeduplicated BAM files are created under the designated directory with *.dedup.bam suffix" >&2
	echo -e "\tIf intermediate or final file already exists, pass the corresponding step or whole steps" >&2
	echo -e "Options:" >&2
        echo -e "\t-o <outDir>: Destination directory, default=<same with the src file>" >&2
        echo -e "\t-m <maxMem>: Maxmum memory, default=5G" >&2
        echo -e "\t-p <thread>: Maximum thread, default=1" >&2
	echo -e "\t-r         : If set, remove duplicate. In default, duplicates are just makred only by 0x400 SAM flag" >&2
#        echo -e "\t-f <samFlag>: SAM file filtering flag, default=NULL" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desBase=NULL
maxMem=5G
thread=1
removeDuplicate=false
while getopts ":o:m:p:r" opt; do
	case $opt in
		o)
			desBase=$OPTARG
			;;
		m)
			maxMem=$OPTARG
			;;
		p)
			thread=$OPTARG
			;;
		r)
			removeDuplicate=true
			;;
#		f)
#			samFlag=$OPTARG
#			;;
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

assertFileExist $@


###################################
## main code

if [ "$desBase" != "NULL" ];then
	mkdir -p $desBase
fi

#maxRead=20000000
memOpt=-Xmx${maxMem}

for src in $@
do
	srcFile=`basename $src`
	if [ "$desBase" == "NULL" ];then
		desDir=`dirname $src`
	else
		desDir=$desBase
	fi

	sorted=${desDir}/${srcFile%.bam}.csort.bam
	dedup=${desDir}/${srcFile%.bam}.dedup.tmp.bam
	final=${desDir}/${srcFile%.bam}.dedup.bam
	metric=${desDir}/${srcFile%.bam}.dedup.metric

	tmp=__temp__.$$.bam
	echo -e "Deduplicate by Picard" >&2
	echo -e "  Src = $src" >&2
	echo -e "  DesSorted = $sorted" >&2
	echo -e "  DesDedup = $final" >&2
	echo -e "  Remove duplicate = ${removeDuplicate}" >&2

	if [ -f $final ];then
		echo -e "  $final already exist, pass" >&2
		echo -e "" >&2
		continue
	fi

	echo -e "  1) Sort by Coordinate => $sorted" >&2
	if [ -f $sorted ];then
		echo -e "  $sorted already exist, pass" >&2
	else
		java $memOpt -jar ${PICARD_HOME}/picard.jar SortSam \
			INPUT=${src} \
			OUTPUT=${tmp} \
			SORT_ORDER=coordinate 
		mv $tmp $sorted
		echo -e "\tSorting $sorted" >&2
		samtools index -b ${sorted} ${sorted}.bai
		#MAX_RECORDS_IN_RAM=$maxRead
	fi

	echo -e "  2) Deduplicate => $dedup" >&2
	if [ -f $dedup ];then
		echo -e "  $dedup already exist, pass" >&2
	else
		java $memOpt -jar ${PICARD_HOME}/picard.jar MarkDuplicates \
			INPUT=${sorted} \
			OUTPUT=${tmp} \
			METRICS_FILE=${metric} \
			REMOVE_DUPLICATES=${removeDuplicate} \
			ASSUME_SORTED=true \
			COMPRESSION_LEVEL=1 
		mv $tmp $dedup
		# MAX_RECORDS_IN_RAM=$maxRead
	fi

	echo -e "  3) Sort by Queryname => $final" >&2
	if [ -f $final ];then
		echo -e "  $final already exist, pass" >&2
	else
		#java $memOpt -jar ${PICARD_HOME}/picard.jar SortSam INPUT=${dedup} OUTPUT=${tmp} SORT_ORDER=queryname MAX_RECORDS_IN_RAM=$maxRead
		samtools sort -n -l 5 -m $maxMem -T __temp__.$$ -@ $thread -o $tmp $dedup
		mv $tmp $final
		rm -v $dedup
	fi

	echo -e "" >&2

done
