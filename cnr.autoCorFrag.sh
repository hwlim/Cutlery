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
	Coordinated-sorted fragment bed file. *.bed.gz or *.bed
Output:
	Three column text file of Distance / Frequency / NormFreq
	NormFreq (Normalized frequency) = Frequency * 1000 / # of total fragment
Options:
	-o <outFile>: Destination directory. Print to stdout if /dev/stdout or stdout. default=stdout
	-m <maxDistance>: maximum distance. default=2000 (bp)
	-v : verbose mode" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=stdout
maxDist=2000
verbose=0
while getopts ":o:m:v" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		m)
			maxDist=$OPTARG
			;;
		v)
			verbose=1
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


if [ $verbose -eq 1 ];then
	echo -e "Calculating fragment autocorrelation" >&2
	echo -e "  src: $src" >&2
	echo -e "  maxDist: $maxDist" >&2
	echo -e "  output: $des" >&2
fi


if [ "$des" == "/dev/stdout" ] || [ "$des" == "stdout" ];then
	des="/dev/stdout"
else
	desDir=`dirname $des`
	mkdir -p $desDir
fi


printBed $src \
	| gawk 'BEGIN{
			maxDist='$maxDist'
			for(i=0;i<=maxDist;i=i+1) dL[i]=0
			pChr=""
			pStart=-1
			ttc=0
		}
		{
			if($1==pChr){
				d=$2-pStart
			}else{
				d=maxDist
				pChr=$1
			}
			pStart=$2
			dL[d]=dL[d]+1
			ttc=ttc+1
		}
		END{
			printf "Distance\tFrequency\tNormFreq\n"
			for(i=0;i<=maxDist;i=i+1) printf "%d\t%d\t%.5f\n", i, dL[i], dL[i]*1000/ttc 
		}' \
	> $des
