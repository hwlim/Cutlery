#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Bam file filtering tool
# - Input: bam files or directories containing multiple bamfiles
#	1. If multiple bam files are given as inputs, individual bam files are processed separately as a distinct sample
#	2. If directories are given as inputs, each directory are treated as single sample
#	   and all bam files within the directory will be merged for filtering
#	Note: BAM file should be sorted by read name (not by coordinate), assuming paired-end sequencing
# - Output: Filtered bam file depending on options,
#	In default, following criteria applied
#	1. Only concordantly aligned pairs (0x2 sam flag)
#	2. Read-pairs with duplicate-flag in sam flag will be discarded (0x400 sam flag)
#	   (no explicit deduplcation is performed in this command)
#	3. Subset of chromosome can be selected such as excluding chrM
#	

source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam file or directory] ..." >&2
	echo -e "Description: Eliminate unnecessary alignment such as chrM or discordant read-pairs" >&2
	echo -e "  1) If a bam file is given as an input:" >&2
	echo -e "\tEach bam file is processed separately" >&2
	echo -e "  2) If a directory is given as an input:" >&2
	echo -e "\tAll bam files under the directory are merged & processed, and single output file is created with the directory name" >&2
	echo -e "Options:" >&2
        echo -e "  -o <outDir>: Destination directory, default=<source directory>" >&2
	echo -e "  -f <samFlag>: SAM flag to include, none if NULL. default=0x2 (properly paired only)" >&2
	echo -e "  -F <samFlag>: SAM flag to exclude, none if NULL. default=0x400 (duplicated)" >&2
	echo -e "  -q <MAPQ>: Alignment quality (MAPQ) threshold for filtering (>= mapq). default=0 (No filtering by MAPQ)" >&2
        echo -e "  -c <chromosome regex>: Regular expression for chromosome selection, default=chr[0-9XY]*$" >&2
        echo -e "                         For multiple patterns use regular expression, such as \"chr[0-9XY]*$|chrM\"" >&2  
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desBase=NULL
flagInc=0x2
flagExc=0x400
chrRegex='chr[0-9XY]*$'
mapq=0
while getopts ":o:f:F:q:c:" opt; do
	case $opt in
		o)
			desBase=$OPTARG
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


printBam(){
	local bam=$1
	assertFileExist $bam

	samtools view -h $optStr $bam
}

printBamDir(){
	local srcDir=$1
	isDirExist ${srcDir}

	srcBamL=( ${srcDir}/*bam )

	samtools view -H ${srcBamL[0]}
	for bam in ${srcBamL[@]}
	do
		samtools view ${optStr} $bam
	done
}

if [ "$desBase" != "NULL" ];then
	mkdir -p $desBase
fi

echo -e "Filtering ATAC-seq data" >&2
echo -e "  SAM flag = $optStr" >&2
echo -e "  MAPQ    >= $mapq" >&2
echo -e "  chrRegex = $chrRegex" >&2
echo -e "  desBase  = $desBase" >&2
echo -e "" >&2

for src in $@
do
	echo -e "Src = $src" >&2
	if [ -f $src ];then
		echo -e "  $src is a file" >&2
		cmd=printBam
		srcName=`basename $src`
		srcName=${srcName%.bam}
	elif [ -d $src ];then
		ls ${src}/*.bam | cut -f 1 | gawk '{ printf "\t - %s\n", $1 }' >&2
		#echo -e "  $src is a directory,\n\tall ${src}/*.bam will be merged" >&2
		cmd=printBamDir
		srcName=`basename $src`
	else
		echo -e "  Error: $src is not a file nor directory" >&2
		exit 1
	fi


	if [ "$desBase" = "NULL" ];then
		desDir=`dirname $src`
	else
		desDir=$desBase
	fi

	des=${desDir}/${srcName}.filtered.bam
	log=${desDir}/${srcName}.filtered.log
	echo -e "Des = $des" >&2

	echo -e "Src = $src" > $log
	echo -e "SAM flag = $optStr" >> $log
	echo -e "chrRegex = $chrRegex" >> $log

	$cmd ${src} \
		| gawk '$3 ~ /'$chrRegex'/ || $1 ~/^@/' \
		| samtools view -b -o __temp__.$$.bam
	mv __temp__.$$.bam $des

	echo -e "" >&2
done
