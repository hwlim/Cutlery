#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bed]
Description:
	Make a raw bedGraph file from a fragment BED file
	**Note that this only considers chromosome names starting with 'chr'
Options:
	-o <outFile>: Destination file. required
	-g <chromSize>: chromosome size file, required
	-m <memory>: memory size for sorting bedGraph file, default=5G
	-c <chrRegex>: regular expression for chromosome filtering, default=. (no filtering)
		e.g. \"^chr[0-9XY]+$|^dm-$\": autosome + sex chromosome + drosophila spikein" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

###################################
## option and input file handling
des=NULL
genome=NULL
memory=5G
chrRegex=.
while getopts ":o:l:L:r:g:m:s:c:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		m)
			memory=$OPTARG
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

assertFileExist $genome

###################################
## main code
printBed(){
	local src=$1
	assertFileExist $src

	local ext=${src##*.}
	
	# avoid using zero coordinate reads
	if [ "$ext" == "gz" ];then
		zcat $src | gawk '$2 > 0'
	else
		cat $src | gawk '$2 > 0'
	fi
}

printFrag(){
	local src=$1
	assertFileExist $src

	printBed $src | gawk '$1 ~ /'$chrRegex'/ { print }'
}

desDir=`dirname $des`
mkdir -p $desDir

tmpBG=${TMPDIR}/__temp__.$$.bedGraph

echo -e "Creating raw bedGraph file from a fragment bed file
  - src = $src
  - des = $des
  - chromSize = $genome
  - chrRegex = $chrRegex" >&2

echo -e "  2) Making bedGraph file" >&2
printFrag $src \
	| sort -S $memory -k1,1 -k2,2n -k3,3n \
	| genomeCoverageBed -bg -g $genome -i stdin \
	| gawk '{ printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4 }' \
	> ${tmpBG}

gzip -c ${tmpBG} > ${des}
rm ${tmpBG}