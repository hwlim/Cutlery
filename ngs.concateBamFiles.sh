#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam1] [bam2] ...
Description:
	Concatenate bam files. No sorting
Options:
	-o <out>: Output file prefix including pathr. required" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=""
while getopts ":o:" opt; do
	case $opt in
		o)
			des=$OPTARG
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


if [ $# -lt 2 ];then
	echo -e "Error: Requires at least two bam files as input" >&2
	exit 1
fi
assertFileExist ${srcL[@]}


if [ "$des" == "" ];then
	echo -e "Error: destination file (-o) must be specified" >&2
	exit 1
fi
desDir=`dirname $des`
mkdir -p $desDir


echo -e "#######################################" >&2
echo -e "Concatenating bam files" >&2
echo -e "  des = $des" >&2
echo -e "  src =" >&2
for src in ${srcL[@]}
do
	echo -e "\t- ${src}" >&2
done

samtools view -H ${srcL[0]} > ${TMPDIR}/__temp__.$$.hdr
samtools cat -h ${TMPDIR}/__temp__.$$.hdr ${srcL[@]} > ${TMPDIR}/__temp__.$$.bam
mv ${TMPDIR}/__temp__.$$.bam ${des}
