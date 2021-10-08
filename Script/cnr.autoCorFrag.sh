#!/usr/bin/env bash

########################################################
## Calculate auto-correlation using CUT&RUN fragments
## - normalized by read depth
##
## To do:
## - homer is very slow, how to make it faster


source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm -rf "${TMPDIR}/__temp__.$$.*"; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bed]
Description: Calculate read fragment auto-correlation by distribution of distance to the next closest fragments
Input:
	Coordinated-sorted fragment bed file. *.bed.gz or *.bed
Output:
	Three column text file of Distance / Frequency / NormFreq
	NormFreq (Normalized frequency) 
		- closest mode: Frequency * 1000000 / # of total fragment
		- homer mode: Frequency * (1000000 / # of total fragment)^2
Options:
	-o <outFile>: Destination directory. Print to stdout if /dev/stdout or stdout. default=stdout
	-m <mode>: auto-correlation mode. 'homer' or 'closest'. default=closest
		in homer mode, all the possible pairs within the maxDistance are counted like homer
		in closest mode, only single closest distance is counted
	-v : verbose mode" >&2
#	-g <chrom>: chromosome size file. chromosomes must be sorted in the same order with the input bed file
#	-d <maxDistance>: maximum distance. default=1000 (bp)
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=stdout
#maxDist=1000
#chrom=""
verbose=0
mode="closest"
while getopts ":o:m:v" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
#		d)
#			maxDist=$OPTARG
#			;;
		m)
			mode=$OPTARG
			;;
#		g)
#			chrom=$OPTARG
#			;;
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

if [ "$mode" != "homer" ] && [ "$mode" != "closest" ];then
	echo -e "Error: mode (-m) must be either homer or closest" >&2
	exit 1
fi

#if [ "$mode" == "homer" ];then
#	if [ "$chrom" == "" ];then
#		echo -e "Error: homer mode requires chrom.size file (-g)" >&2
#		exit 1
#	else
#		assertFileExist $chrom
#	fi
#fi

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

maxDist=1999
if [ $verbose -eq 1 ];then
	echo -e "Calculating fragment autocorrelation" >&2
	echo -e "  src: $src" >&2
	echo -e "  mode: $mode" >&2
	echo -e "  maxDist: $maxDist" >&2
	echo -e "  output: $des" >&2
fi


if [ "$des" == "/dev/stdout" ] || [ "$des" == "stdout" ];then
	des="/dev/stdout"
else
	desDir=`dirname $des`
	mkdir -p $desDir
fi

if [ "$mode" == "closest" ];then
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
else
	ttc=`printBed $src | wc -l`
	tmpDir=${TMPDIR}/__temp__.$$.homer
	printBed $src \
		| makeTagDirectory $tmpDir /dev/stdin 2>/dev/null
	
	tail -n +2 ${tmpDir}/tagAutocorrelation.txt \
		| gawk 'BEGIN{
					normFactor = (1000000/'$ttc') ^ 2
					printf "Distance\tFrequency\tNormFreq\n"
				}
				{
					if($1 >= 0) printf "%d\t%d\t%.3f\n", $1,$2,$2*normFactor
				}' \
		> $des
	rm -rf $tmpDir
	#printBed $src \
	#	| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1,$2,$2+'$maxDist',$4,$5,$6 }' \
	#      	| intersectBed -a stdin -b $src -wa -wb -sorted \
	#	| gawk 'BEGIN{
	#			ttc='$ttc'
	#			maxDist='$maxDist'
	#			for(i=0;i<=maxDist;i=i+1) dL[i]=0
	#		}
	#		{
	#			if($4==$10) next
	#			d = $8 - $2
	#			dL[d] = dL[d] + 1
	#		}
	#		END{
	#			printf "Distance\tFrequency\tNormFreq\n"
	#			for(i=0;i<=maxDist;i=i+1) printf "%d\t%d\t%.5f\n", i, dL[i], dL[i]*1000/ttc
	#		}' \
	#	> $des
fi
