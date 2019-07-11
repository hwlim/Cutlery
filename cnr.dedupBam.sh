#!/usr/bin/env bash

# ATAC-seq tools
# Written by Hee-Wooong Lim
# 
source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam]" >&2
	echo -e "Description:" >&2
	echo -e "\tDeduplicate given bam files using Picard tool" >&2
	echo -e "\t1) Coordinate sorting => 2) Deduplication => 3) Readname sorting" >&2
	echo -e "\tIf intermediate or final file already exists, pass the corresponding step or whole steps" >&2
	echo -e "Options:" >&2
	echo -e "\t-o <outFile>: output files, default=<src file name without path/.bam>.dedup.bam" >&2
        echo -e "\t-m <maxMem>: Maxmum memory, default=5G" >&2
#        echo -e "\t-p <thread>: Maximum thread, default=1" >&2
	echo -e "\t-r         : If set, remove duplicate. In default, duplicates are just makred only by 0x400 SAM flag" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=NULL
maxMem=5G
#thread=1
removeDuplicate=false
while getopts ":o:m:r" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		m)
			maxMem=$OPTARG
			;;
#		p)
#			thread=$OPTARG
#			;;
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

src=$1
assertFileExist $src


###################################
## main code


#maxRead=20000000
memOpt=-Xmx${maxMem}

if [ "$des" = "NULL" ];then
	srcPrefix=`basename $src .bam`
	des=${srcPrefix}.dedup.bam
fi
desDir=`dirname $des`
mkdir -p $desDir

## coordindate sorted bam
## deduplicated / coordinate sorted bam
## deduplicated / readname sorted bam
outPrefix=${des%.bam}
tmpCsort=${TMPDIR}/__temp__.$$.csort.bam
tmpDedupCsort=${TMPDIR}/__temp__.$$.csort.dedup.bam
tmpDedup=${TMPDIR}/__temp__.$$.dedup.bam
metric=${outPrefix}.metric

if [ "$src" == "$des" ];then
	echo -e "Error: source and destination are the same" >&2
	exit 1
fi
tmp=${TMPDIR}/__temp__.$$.bam
echo -e "Deduplicate by Picard" >&2
echo -e "- Src = $src" >&2
echo -e "- Des = $des" >&2
echo -e "- Remove duplicate = ${removeDuplicate}" >&2

echo -e "  1) Sort by Coordinate => $tmpCsort" >&2
java $memOpt -jar ${PICARD_HOME}/picard.jar SortSam \
	INPUT=${src} \
	OUTPUT=${tmpCsort} \
	SORT_ORDER=coordinate \
	COMPRESSION_LEVEL=1 

echo -e "  2) Deduplicate => $tmpDedupCsort" >&2
java $memOpt -jar ${PICARD_HOME}/picard.jar MarkDuplicates \
	INPUT=${tmpCsort} \
	OUTPUT=${tmpDedupCsort} \
	METRICS_FILE=${metric} \
	REMOVE_DUPLICATES=${removeDuplicate} \
	ASSUME_SORTED=true \
	COMPRESSION_LEVEL=1 

echo -e "  3) Sort by Queryname => $des" >&2
java $memOpt -jar ${PICARD_HOME}/picard.jar SortSam \
	INPUT=${tmpDedupCsort} \
	OUTPUT=${tmpDedup} \
	SORT_ORDER=queryname 
#java $memOpt -jar ${PICARD_HOME}/picard.jar SortSam INPUT=${dedup} OUTPUT=${tmp} SORT_ORDER=queryname MAX_RECORDS_IN_RAM=$maxRead
#samtools sort -n -l 5 -m $maxMem -T ${TMPDIR}/__temp__.$$ -@ $thread -o $tmpDedup $tmpDedupCsort
mv $tmpDedup $des

