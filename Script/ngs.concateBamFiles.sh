#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam1] ...
Description:
	Merge multiple bam files into one
Input:
	- Sorted or unsorted BAM files
Output:
	- Simply concatenated BAM file without indexing (default)
	or
	- Coordinate-sorted bam file with index, <out>.bai (with -S or -s option)
Options:
	-o <out>: Output bam file including path. required.
		Index will be automatically named as <out>.bai if created
	-m <mem>: Memory for sorting (for -s option). default=5G
	-S : If set, input bam files are assumed to be coordinated-sorted (index files must exist)
		So, input files are merged by \"samtools merge\".
		Hence the result is also coordinate-sorted and indexed
	-s : If set, perform coordinate-sorting and indexing assuming unsorted input bam files.
		If -S is set, this option is ignored" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=""
sortMem=5G
sortedBam=FALSE
toSort=FALSE
while getopts ":o:m:Ss" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		m)
			sortMem=$OPTARG
			;;
		S)
			sortedBam=TRUE
			;;
		s)
			toSort=TRUE
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
echo -e "  sorted input = $sortedBam" >&2
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
	if [ "$sortedBam" == "TRUE" ];then
		samtools index $des
	fi
else
	if [ "$sortedBam" == "TRUE" ];then
		echo -e "- Input bam files are pre-sorted. Performing samtools merge & indexing" >&2
		samtools merge $tmpMerged ${srcL[@]}
		samtools index $tmpMerged
		mv $tmpMerged $des
		mv ${tmpMerged}.bai ${des}.bai
	else
		echo -e "- Simply concatenating input bam files" >&2
		samtools view -H ${srcL[0]} > ${tmpHdr}
		samtools cat -h ${tmpHdr} ${srcL[@]} > ${tmpMerged}

		if [ "$toSort" == "TRUE" ];then
			echo -e "- Sorting & Indexing" >&2
			samtools sort -o $tmpSrt -T $tmpSrtPrefix -m $sortMem $tmpMerged
			mv $tmpSrt $des
			samtools index $des
		else
			mv $tmpMerged $des
		fi
	fi
fi

