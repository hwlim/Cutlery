#!/usr/bin/env bash


##############################################################################
# Written by Hee-Wooong Lim
# - Convert given alignment file (bam or fragment bed) into v-plot coordiate:
#	chr / start / end / name / fragLen / +
#	where start = end = fragment center
	
source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [bam/bed.gz]
Description: Convert paired-end BAM file or fragment bed file into a bed file containing center/fragLen
	BAM file: paired-end & name-sorted file 
	BED file: gzipped bed file where each line correspond to single fragment
Options:
	-o <outFile>: output file. Sorted by coordinate. default=<src file with path>.cfl.bed.gz
		chr / start / end / name / fragLen / +
		where start = end = fragment center
	-m: Memory size for sort, e.g. 10G. default=5G" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=NULL
sortMem=5G
while getopts ":o:m:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		m)
			sortMem=$OPTARG
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

src=$1
assertFileExist $src


###################################
## main code

if [ "$des" == "NULL" ];then
	srcFile=`basename $src`
	des=${srcFile%.bam}.frag.bed.gz
fi
if [ "$des" == "$src" ];then
	echo -e" Error: source and destination files are the same" >&2
	exi 1
fi

desDir=`dirname $des`
mkdir -p $desDir


tmpDes=${TMPDIR}/__temp__.$$.bed.gz

#################################
## Processing
ext=${src##*.}
printFrag(){
	if [ "$ext" == "gz" ];then
		zcat $src
	elif [ "$ext" == "bam" ];then
		bamToBed -bedpe -i $src \
			| gawk '{ printf "%s\t%d\t%d\tFrag.%d\t0\t+\n", $1,$2,$6,NR }'
	fi
}


echo -e "=======================" >&2
echo -e "BAM to fragment" >&2
echo -e "- src = $src" >&2
echo -e "- des = $des" >&2
echo -e "- sortMem = $sortMem" >&2

printFrag $src \
	| gawk '{
			c=($2+$3)
			printf "%s\t%d\t%d\t%s\t%d\t+\n", $1,c,c,$4,$3-$2 
		}' \
	| sort -S $sortMem -k1,1 -k2,2n -k3,3n \
	| gzip \
	> $tmpDes

mv $tmpDes ${des}
