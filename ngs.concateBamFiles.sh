#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam1] ...
Description:
	Merge multiple bam files into one
	In default, bam files are simply concatenated
Options:
	-o <out>: Output file prefix including path. required
	-m <mem>: Memory for sorting. default=5G
	-s : If set, coordinate-sorted and indexed" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=""
sortMem=5G
toSort=0
while getopts ":o:m:s" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		m)
			sortMem=$OPTARG
			;;
		s)
			toSort=1
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
srcL=( $@ )


if [ $# -lt 1 ];then
	echo -e "Error: Requires at least one bam file as input" >&2
	exit 1
fi

if [ "$des" == "" ];then
	echo -e "Error: destination file (-o) must be specified" >&2
	exit 1
fi

assertFileExist ${srcL[@]}
desDir=`dirname $des`
mkdir -p $desDir

echo -e "#######################################" >&2
echo -e "Concatenating bam files" >&2
echo -e "  des = $des" >&2
echo -e "  sort = $toSort" >&2
echo -e "  memory = $sortMem" >&2
echo -e "  src =" >&2
for src in ${srcL[@]}
do
	echo -e "\t- ${src}" >&2
done

tmpHdr=${TMPDIR}/__temp__.$$.hdr
tmpMerged=${TMPDIR}/__temp__.$$.bam
tmpSrt=${TMPDIR}/__temp__.$$.sorted.bam
tmpSrtPrefix=${TMPDIR}/__temp__.$$.sorted

if [ $# -eq 1 ];then
	echo -e "Warning: Only one file is given; simply copying" >&2
	cp ${srcL[0]} $des
else
	samtools view -H ${srcL[0]} > ${tmpHdr}
	samtools cat -h ${tmpHdr} ${srcL[@]} > ${tmpMerged}

	if [ $toSort -gt 0 ];then
		samtools sort -o $tmpSrt -T $tmpSrtPrefix -m $sortMem $tmpMerged
		mv $tmpSrt $des
		samtools index $des
	else
		mv $tmpMerged $des
	fi
fi

