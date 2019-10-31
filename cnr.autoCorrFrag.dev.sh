#!/usr/bin/env bash

########################################################
## Calculate auto-correlation using CUT&RUN fragments
## - normalized by read depth


source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bed]
Description: Calculate read fragment auto-correlation by distribution of distance to the next closest fragments
Input:
	Coordinated-sorted fragment bed file
Output:
	Three column text file of Distance / Frequency / NormalizedFreq
	Normalized frequency = Frequency * 1000 / total fragment
Options:
	-o <outFile>: Destination directory. Print to stdout if /dev/stdout or stdout. default=stdout
	-m <maxDistance>: maximum distance. default=2000 (bp)" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=stdout
maxDist=2000
while getopts ":o:m:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		m)
			memory=$OPTARG
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

if [ "$genome" = "NULL" ];then
	echo -e "Error: genome (-g) must be specified" >&2
	exit 1
fi
src=$1
assertFileExist $src

if [ "$des" == "NULL" ];then
	echo -e "Error: Destination file must be specified" >&2
	exit 1
fi

if [ -f $genome ];then
	chrom=${genome}
else
	chrom=~/Research/Common_Data/${genome}/chrom.sizes
fi

###################################
## main code
printBed(){
	local src=$1
	assertFileExist $src

	local ext=${src##*.}
	
	if [ "$ext" == "gz" ];then
		zcat $src | grep ^chr
	else
		cat $src | grep ^chr
	fi
}

assertFileExist $chrom

desDir=`dirname $des`
mkdir -p $desDir

tmpBG=${TMPDIR}/__temp__.$$.bedGraph
tmpBW=${TMPDIR}/__temp__.$$.bw

echo -e "Creating BigWig file from a fragment bed file" >&2
echo -e "- src = $src" >&2
echo -e "- des = $des" >&2
echo -e "- genome = $genome" >&2

if [ $scaleFactor == "0" ];then
	echo -e "  1) Calculating scale factor for RPM normalization" >&2
	ttc=`printBed $src | wc -l`
	scaleFactor=`echo $ttc | gawk '{ printf "%f", 1000000/$1}'`
	echo -e "\tTTC = $ttc (scaleFactor $scaleFactor)" >&2
else
	echo -e "  1) Scale factor was manually assigned" >&2
	echo -e "\tscaleFactor: ${scaleFactor}" >&2
fi

echo -e "  2) Making bedGraph file" >&2
printBed $src \
	| sort -S $memory -k1,1 -k2,2n -k3,3n \
	| genomeCoverageBed -bg -scale $scaleFactor -g $chrom -i stdin \
	| gawk '{ printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4 }' \
	> $tmpBG

echo -e "  3) Converting to bigWig file" >&2
bedGraphToBigWig ${tmpBG} $chrom ${tmpBW}
mv ${tmpBW} ${des}

