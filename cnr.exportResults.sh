#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <src dir> <des dir>
Description: Export Cutlery processing results to a destination folder
Input:
	A folder containing sample folders that containing outputs, i.e. stratified by samples
	e.g.
	3.Sample/
	├── Ctrl
	├── dnRAR
	└── v2
Output:
	A folder containing outputs organized by output types and prefixed by sample name
Options:
	-f : Overwrite existing files/folders if set" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi



###################################
## option and input file handling
force=FALSE
#lengthParam=NULL,NULL
while getopts ":f" opt; do
	case $opt in
		f)
			force=TRUE
			;;
#		s)
#			optStr=$OPTARG
#			;;
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
if [ $# -lt 2 ];then
	printUsage
	exit 1
fi

srcDir=$1
desDir=$2

assertDirExist $srcDir
mkdir -p $desDir


echo -e "Exporting Cutlery results" >&2
echo -e "  src: $srcDir" >&2
echo -e "  des: $desDir" >&2
echo -e "  force: $force" >&2


sampleL=( `ls -d ${srcDir}/*/` )

## export list
# - bigwig: all / nfr / nuc
# - peak & heatmap: all / 1rpm
# - motif

exportFile(){
	local src=$1
	local des=$2
	local isDir=$3

	if [ "$isDir" == "TRUE" ];then
		optStr="-r"
	fi

	if [ -f $src ];then
		if [ -f $des ];then
		       if [ "$force" == "TRUE" ];then
			       echo -e "Warning: Overwriting $des" >&2
				cp $optStr -f -v $src $des
			else
				echo -e "Warning: $des already exists, skip" >&2
			fi
		else
			cp $optStr -v $src $des
		fi
	else
		echo -e "Warning: $src does not exist, skipping" >&2
	fi
}

for samplePath in ${sampleL[@]}
do
	sample=`basename $samplePath`

	# fragment length distribution
	echo -e "0) Exporting fragment length distribution" >&2
	mkdir -p ${desDir}/QualityControl
	src=${srcDir}/${sample}/QC/fragLen.dist.png
	des=${desDir}/QualityControl/${sample}.fragLen.png
	exportFile $src $des FALSE
	#if [ -f $src ];then
	#	if [ ! -f $des ] || [ "$force" == "TRUE" ];then
	#		cp -f -v $src $des
	#	else
	#		echo -e "Warning: $des already exists, ski" >&2
	#	fi
	#else
	#	echo -e "Warning: $src does not exist, skipping" >&2
	#fi


	# - bigwig: all / nfr / nuc
	echo -e "1) Exporting bigwig files" >&2
	mkdir -p ${desDir}/BigWig
	for frag in all nfr nuc
	do
		src=${srcDir}/${sample}/igv.${frag}.ctr.bw
		des=${desDir}/BigWig/${sample}.${frag}.ctr.bw
		exportFile $src $des FALSE

		#if [ -f $src ];then
		#	if [ ! -f $des ] || [ "$force" == "TRUE" ];then
		#		cp -f -v $src $des
		#	else
		#		echo -e "Warning: $des already exists, ski" >&2
		#	fi
		#else
		#	echo -e "Warning: $src does not exist, skipping" >&2
		#fi
	done

	# - peak & heatmap: all / 1rpm
	echo -e "2) Exporting peak files & heatmaps" >&2
	mkdir -p ${desDir}/Peak
	for mode in factor histone
	do
		src=${srcDir}/${sample}/HomerPeak.${mode}/peak.exBL.bed
		des=${desDir}/Peak/${sample}.all.bed
		exportFile $src $des FALSE
		#if [ -f $src ];then
		#	if [ ! -f $des ] || [ "$force" == "TRUE" ];then
		#		cp -f -v $src $des
		#	else
		#		echo -e "Warning: $des already exists, ski" >&2
		#	fi
		#else
		#	echo -e "Warning: $src does not exist, skipping" >&2
		#fi

		src=${srcDir}/${sample}/HomerPeak.${mode}/peak.exBL.1rpm.bed
		des=${desDir}/Peak/${sample}.1rpm.bed
		exportFile $src $des FALSE
		#if [ -f $src ];then
		#	if [ ! -f $des ] || [ "$force" == "TRUE" ];then
		#		cp -f -v $src $des
		#	else
		#		echo -e "Warning: $des already exists, ski" >&2
		#	fi
		#else
		#	echo -e "Warning: $src does not exist, skipping" >&2
		#fi

		src=${srcDir}/${sample}/HomerPeak.${mode}/heatmap.exBL.1rpm.png
		des=${desDir}/Peak/${sample}.1rpm.heatmap.png
		exportFile $src $des FALSE
		#if [ -f $src ];then
		#	if [ ! -f $des ] || [ "$force" == "TRUE" ];then
		#		cp -f -v $src $des
		#	else
		#		echo -e "Warning: $des already exists, ski" >&2
		#	fi
		#else
		#	echo -e "Warning: $src does not exist, skipping" >&2
		#fi
	done

	# - motif
	echo -e "3) Exporting motif search results" >&2
	mkdir -p ${desDir}/Motif
	for src in ${srcDir}/${sample}/Motif/Homer.all  ${srcDir}/${sample}/Motif/MEME.random5k
	do
		mode=`basename $src`
		des=${desDir}/Motif/${sample}.${mode}
		exportFile $src $des TRUE
		#if [ -d $src ];then
		#	if [ ! -d $des ] || [ "$force" == "TRUE" ];then
		#		cp -f -rv $src $des
		#	else
		#		echo -e "Warning: $des already exists, ski" >&2
		#	fi
		#else
		#	echo -e "Warning: $src does not exist, skipping" >&2
		#fi
	done
done
