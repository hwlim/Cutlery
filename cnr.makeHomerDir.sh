#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 
# Wrapper script to make Homer tag directory
#
# - Input: fragment BED file (*nfr.ct.bed) 
# - Output: 
#	<oudesDir>/TSV	: Homer tag directory
#	<oudesDir>/HomerPeak : peak directory
#

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [fragment bed file in gzip]
Description: Make Homer data directory from BED file
Options:
	-o <outDir>: Destination tag directory, required
	-n <name>: data name to be stored under */TSV/info.txt, default=<src file name>
	-c <chrRegex>: regular expression for chromosome filtering, default=NULL (no filtering)
		\"^chr[0-9XY]+$|^dm-$\": autosome + sex chromosome + drosophila spikein" >&2
#	echo -e "\t-l <fragLen,peakWidth>: Comma-separated fragment length and peak width, default=NULL,NULL" >&2
#        echo -e "\t-g <genome>: genome, default=NULL" >&2
#        echo -e "\t-h: Print help" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desDir=NULL
name=NULL
chrRegex=NULL
#lengthParam=NULL,NULL
while getopts ":o:n:c:" opt; do
	case $opt in
		o)
			desDir=$OPTARG
			;;
		n)
			name=$OPTARG
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
assertFileExist $src
if [ "$desDir" == "NULL" ];then
	echo -e "Error: Destination directory (-o) must be specified" >&2
	exit 1
fi

if [ "$name" == "NULL" ];then
	name=`basename $src`
fi



###################################
## main code
printFile(){
	local src=$1
	assertFileExist $src

	local ext=${src##*.}
	
	if [ "$ext" == "gz" ];then
		zcat $src
	else
		cat $src
	fi
}



tmpTagDir=${TMPDIR}/__temp__.$$.tDir
log=${tmpTagDir}/TSV.log
echo -e "Creating Homer tag directory" >&2
echo -e "- src = $src" >&2
echo -e "- name = $name" >&2
echo -e "- desDir = $desDir" >&2
echo -e "- chrRegex = $chrRegex" >&2

#mkdir -p $desDir
echo -e "$name" > ${tmpTagDir}/info.txt

if [ "$chrRegex" == "NULL" ];then
	printFile $src \
		| makeTagDirectory ${tmpTagDir} /dev/stdin -format bed -fragLength given 2>&1 | tee $log
else
	printFile $src \
		| gawk '$1 ~ /'$chrRegex'/' \
		| makeTagDirectory ${tmpTagDir} /dev/stdin -format bed -fragLength given 2>&1 | tee $log
fi

drawAutoCorrplot.r -t "$name" ${tmpTagDir}
mv ${tmpTagDir} ${desDir}

