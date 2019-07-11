#!/usr/bin/env bash

# command-line application template

source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bed]" >&2
	echo -e "Description: Make a bigWig file from fragment BED file" >&2
	echo -e "Options:" >&2
        echo -e "\t-o <outFile>: Destination directory. required" >&2
        echo -e "\t-g <genome>: genome or chromosome size file, default=NULL" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=NULL
genome=NULL
while getopts ":o:g:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		g)
			genome=$OPTARG
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
desDir=`dirname $des`
mkdir -p $desDir


###################################
## main code
printBed(){
	local src=$1
	assertFileExist $src

	local ext=${src##*.}
	
	if [ "$ext" == "gz" ];then
		zcat $src
	else
		cat $src
	fi
}

if [ -f $genome ];then
	chrom=${genome}
else
	chrom=~/Research/Common_Data/${genome}/chrom.sizes
fi
assertFileExist $chrom

srcFile=`basename $src`
if [ "$desDir" == "NULL" ];then
	desDir=`dirname $src`
fi
tmpBG=__temp__.$$.bedGraph
tmpBW=__temp__.$$.bw
des=${desDir}/${srcFile%.bed*}.bw

echo -e "Creating BigWig file from a fragment bed file" >&2
echo -e "- src = $src" >&2
echo -e "- des = $des" >&2
echo -e "- genome = $genome" >&2

echo -e "  1) Counting total tag count" >&2
ttc=`printBed $src | wc -l`
scale=`echo $ttc | gawk '{ printf "%f", 1000000/$1}'`
echo -e "\tTTC = $ttc (scale $scale)" >&2

echo -e "  2) Making bedGraph file" >&2
printBed $src \
	| sort -S 10G -k1,1 -k2,2n -k3,3n \
	| genomeCoverageBed -bg -scale $scale -g $chrom -i stdin \
	| gawk '{ printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4 }' \
	> $tmpBG

echo -e "  3) Converting to bigWig file" >&2
bedGraphToBigWig ${tmpBG} $chrom ${tmpBW}
mv ${tmpBW} ${des}

