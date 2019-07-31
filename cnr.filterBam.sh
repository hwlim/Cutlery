#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Bam file filtering tool
# - Input: bam file(s)
# - Output: Filtered bam file depending on options,
#	In default, following criteria applied
#	1. Only concordantly aligned pairs (0x2 sam flag)
#	2. Read-pairs with duplicate-flag in sam flag will be discarded (0x400 sam flag)
#	   (no explicit deduplcation is performed in this command)
#	3. Subset of chromosome can be selected such as excluding chrM
#	
# Technical consideration
# - chromosome selection in samtools view only works for sorted bam file
#   So, for unsorted file, chromosome must be checked explicitly

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam]
Description:
	Eliminate unnecessary alignment such as chrM or discordant read-pairs
	If multiple bam files are given, they are processed all together and output is saved in a single bam file
Input:
	Paired-end BAM file(s)
Options:
        -o <output>: output file, default=<src file name>.filtered.bam
	-f <samFlag>: SAM flag to include, none if NULL. default=0x2 (properly paired only)
	-F <samFlag>: SAM flag to exclude, none if NULL. default=0x400 (duplicated)
	-q <MAPQ>: Alignment quality (MAPQ) threshold for filtering (>= mapq). default=0 (No filtering by MAPQ)
        -c <chromosome regex>: Regular expression for chromosome selection, default=^chr[0-9XY]+$
		For multiple patterns use regular expression, such as \"^chr[0-9XY]+$|chrM\" 
		NULL if not applicable or no filtering" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=""
flagInc=0x2
flagExc=0x400
chrRegex='^chr[0-9XY]+$'
mapq=0
while getopts ":o:f:F:q:c:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		f)
			flagInc=$OPTARG
			;;
		F)
			flagExc=$OPTARG
			;;
		q)
			mapq=$OPTARG
			;;
		c)
			chrRegex=$OPTARG
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


###################################
## main code

optStr=""
if [ "$flagInc" != "NULL" ];then
	optStr="-f $flagInc"
fi
if [ "$flagExc" != "NULL" ];then
	optStr="${optStr} -F $flagExc"
fi
if [ $mapq -gt 0 ];then
	optStr="${optStr} -q $mapq"
fi

assertFileExist $src
srcFile=`basename $src`
if [ "$des" == "" ];then
	des=${desFile%.bam}.filtered.bam
fi
if [ "$des" == "$src" ];then
	echo -e "Error: source and destination files are the same" >&2
	exit 1
fi
desLog=${des%.bam}.log

desDir=`dirname $des`
mkdir -p $desDir


echo -e "Filtering BAM file" >&2
echo -e "- SAM flag = $optStr" >&2
echo -e "- MAPQ    >= $mapq" >&2
echo -e "- chrRegex = $chrRegex" >&2
echo -e "- src  = $src" >&2
echo -e "- des  = $des" >&2

echo -e "Filtering BAM file" > $desLog
echo -e "- SAM flag = $optStr" >> $desLog
echo -e "- MAPQ    >= $mapq" >> $desLog
echo -e "- chrRegex = $chrRegex" >> $desLog
echo -e "- src  = $src" >> $desLog
echo -e "- des  = $des" >> $desLog


tmp=${TMPDIR}/__temp__.$$.bam
#chrList=`samtools view -H $src | sed 's/:/\t/' | gawk '{ if($1=="@SQ" && $2=="SN") print $3 }' | grep -E -w ${chrRegex}`
#samtools view -b -o $tmp $src $chrList 

if [ "$chrRegex" == "NULL" ];then
	samtools view -b -1 optStr -o $tmp $src
else
	samtools view -h $optStr $src \
		| gawk '$3 ~ /'$chrRegex'/ || $1 ~/^@/' \
		| samtools view -b -1 -o $tmp
fi
mv $tmp $des
