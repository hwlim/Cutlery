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
	-l <minLen>: lower bound of fragment length to use (including the boundary), default=0
	-L <maxLen>: upper bound of fragment length to use (including the boundary), default=1000000 (infinite)
	-r <resize>: resize to the give length around the center, default=-1 (no resize)
	-n <name>: data name to be stored in */TSV/info.txt, default=<src file name>
	-c <chrRegex>: regular expression for chromosome filtering, default=. (no filtering)
		e.g. \"^chr[0-9XY]+$|^dm-$\": autosome + sex chromosome + drosophila spikein" >&2
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
minLen=0
maxLen=1000000
resize=-1
name=NULL
chrRegex=NULL
#lengthParam=NULL,NULL
while getopts ":o:l:L:r:n:c:" opt; do
	case $opt in
		o)
			desDir=$OPTARG
			;;
		l)
			minLen=$OPTARG
			;;
		L)
			maxLen=$OPTARG
			;;
		r)
			resize=$OPTARG
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


printFrag(){
	local src=$1
	assertFileExist $src

	printFile $src \
		| gawk 'BEGIN{
				minLen='$minLen';
				maxLen='$maxLen';
				resize='$resize';
			}{
				fragLen=$3-$2
				if( fragLen < minLen || fragLen > maxLen ) next

				if( resize > -1 ){
					c = ($3+$2)/2
					h = resize / 2
					printf "%s\t%d\t%d\tNULL\t0\t%s\n", $1, c-h, c+h, $6
				}else{
					printf "%s\t%d\t%d\tNULL\t0\t%s\n", $1, $2, $3, $6
				}
			}'
}


tmpTagDir=${TMPDIR}/__temp__.$$.tDir
log=${tmpTagDir}/TSV.log
echo -e "Creating Homer tag directory
  - src = $src
  - fragLen = $minLen - $maxLen (bp)
  - resize = $resize (bp)
  - name = $name
  - desDir = $desDir
  - chrRegex = $chrRegex" >&2

mkdir -p $tmpTagDir
echo -e "$name" > ${tmpTagDir}/info.txt

#if [ "$chrRegex" == "NULL" ];then
#	printFrag $src \
#		| makeTagDirectory ${tmpTagDir} /dev/stdin -format bed -fragLength given 2>&1 | tee $log
#else
	printFrag $src \
		| gawk '$1 ~ /'$chrRegex'/' \
		| makeTagDirectory ${tmpTagDir} /dev/stdin -format bed -fragLength given 2>&1 | tee $log
#fi

drawHomerAutoCorr.r -t "$name" ${tmpTagDir}
mv ${tmpTagDir} ${desDir}

