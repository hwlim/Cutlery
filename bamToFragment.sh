#!/usr/bin/env bash


##############################################################################
# Written by Hee-Wooong Lim
# - Convert paired-end BAM file into fragment bed file by connecting two reads
# - Assume that reads are sorted by query name
#
	
source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [bam]
Description: Convert paired-end BAM file into fragment bed file by connecting two reads
	Assume that reads reads are sorted by query name
Options:
	-o <outFile>: output file.  default=<bamFile without path & extension>.frag.bed.gz
		e.g. ../input.bam -> input.frag.bed.gz
	-l <fragLen>: resize the fragment around the center. 0 for no resize. default=0
	-s: Sort by sort -k1,1 -k2,2n -k3,3n. default=Off
	-m: Memory size for sort, e.g. 10G. default=5G" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=NULL
fragLen=0
sortBed=FALSE
sortMem=5G
while getopts ":o:l:m:s" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		l)
			fragLen=$OPTARG
			;;
		m)
			sortMem=$OPTARG
			;;
		s)
			sortBam=TRUE
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


#optStr=""
#if [ "$flagInc" != "NULL" ];then
#	optStr="-f $flagInc"
#fi
#if [ "$flagExc" != "NULL" ];then
#	optStr="$optStr -F $flagExc"
#fi

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


#################################
## Processing

echo -e "=======================" >&2
echo -e "BAM to fragment" >&2
echo -e "- src = $src" >&2
echo -e "- des = $des" >&2
echo -e "- fragLen = $fragLen" >&2
echo -e "- sortMem = $sortMem" >&2
echo -e "- TMPDIR = $TMPDIR" >&2

tmpDes=${TMPDIR}/__temp__.$$.bed.gz


printBed(){
	if [ $fragLen -eq 0 ];then
		bamToBed -bedpe -i $1 \
			| gawk '{ printf "%s\t%d\t%d\tFrag.%d\t0\t+\n", $1,$2,$6,NR }'
	else
		bamToBed -bedpe -i $1 \
			| gawk 'BEGIN{
					hWid='$fragLen'/2
				}
				{
					d=$6-$2
					c=($6+$2)/2
					printf "%s\t%d\t%d\tFrag.%d\t0\t+\n", $1,c-hWid,c+hWid,NR
				}'
	fi
}

if [ "$sortBam" = "TRUE" ];then
	printBed $src \
		| sort -S $sortMem -k1,1 -k2,2n -k3,3n \
		| gzip \
		> $tmpDes
else
	printBed $src \
		| gzip \
		> $tmpDes
fi

mv $tmpDes ${des}
