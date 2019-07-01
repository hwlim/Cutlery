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

source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam1] [bam2] ..." >&2
	echo -e "Description:" >&2
	echo -e "\tEliminate unnecessary alignment such as chrM or discordant read-pairs" >&2
	echo -e "\tIf multiple bam files are given, they are processed all together and output is saved in a single bam file" >&2
	echo -e "Input:" >&2
	echo -e "\tPaired-end BAM file(s)" >&2
	echo -e "Options:" >&2
        echo -e "\t-o <output>: output file, must be specified" >&2
	echo -e "\t-f <samFlag>: SAM flag to include, none if NULL. default=0x2 (properly paired only)" >&2
	echo -e "\t-F <samFlag>: SAM flag to exclude, none if NULL. default=0x400 (duplicated)" >&2
	echo -e "\t-q <MAPQ>: Alignment quality (MAPQ) threshold for filtering (>= mapq). default=0 (No filtering by MAPQ)" >&2
        echo -e "\t-c <chromosome regex>: Regular expression for chromosome selection, default=chr[0-9XY]*$" >&2
        echo -e "\t\tFor multiple patterns use regular expression, such as \"chr[0-9XY]*$|chrM\"" >&2  
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
chrRegex='chr[0-9XY]*$'
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


desDir=`dirname $des`
mkdir -p $desDir

srcL=( $@ )

echo -e "Filtering BAM file(s)" >&2
echo -e "  SAM flag = $optStr" >&2
echo -e "  MAPQ    >= $mapq" >&2
echo -e "  chrRegex = $chrRegex" >&2
echo -e "  des  = $des" >&2

chrList=`samtools view -H ${srcL[0]} | grep ... `

tmpL=""
for (( i=0;i<${#srcL[@]};i=i+1 ))
do
	echo -e "  - Filtering : $src" >&2
	tmp=${TMPDIR}/__temp__.$$.${i}.bam
	tmpL="${tmpL} ${tmp}"

	samtools view -b -o $tmp $chrList > ${tmp}
done

tmpHdr=${TMPDIR}/__temp__.$$.hdr.sam
samtools view -H ${srcL[0]} > $tmpHdr
if [ ${#srcL[@]} -gt 1 ];then
	echo -e "  - Merging" >&2
	tmpDes=${TMPDIR}/__temp__.$$.bam
	samtools merge -h $tmpHdr $tmpDes $tmpL
	mv $tmpDes $des
else
	mv $tmpL $des
fi
