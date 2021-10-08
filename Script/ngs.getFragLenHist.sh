#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [bam or bed.gz]
Description:
	- Check the frequency of fragment length
	- *Warning: Chromosomes start with 'chr' are considered.
	- Fragment length is the distance between 5'-ends of Read1 and Read2, i.e. handling protrusion properly
            |------R1------->
          |<-------R2-----|
            |--- length---|
Input:
	- Paired-end bam file (assuming sorted-by-name in default)
		Coordinate-sorted bam file is also allowed (-c option)
	- Fragment bed file (compressed)
Output:
	Two column text file with header, 'fragLen' / 'Cnt'
Options:
	-o <out>: output file, default=stdout
	-l <maxLen>: max length to consider. Larger lengths are included in the maxLen, default=1000
	-c : If set, input bam file is assumed to be coordinate-sorted"
}


###################################
## option and input file handling
out=NULL
maxLen=1000
cSorted=FALSE
#chregex=
while getopts ":o:l:c" opt; do
	case $opt in
		o)
			out=$OPTARG
			;;
		l)
			maxLen=$OPTARG
			;;
		c)
			cSorted=TRUE
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

if [ "$out" = "NULL" ];then
	out=/dev/stdout
else
	desDir=`dirname $out`
	mkdir -p $desDir
fi

###################################
## main code

ext=${src##*.}

if [ $ext = "bam" ];then
	if [ "$cSorted" = "TRUE" ];then
		sortedBamToFrag.py $src \
			| grep ^chr \
			| gawk 'BEGIN{ printf "fragLen\tCnt\n"; maxLen='$maxLen' }
				{
					d=$3-$2
					if( d>maxLen ){ d=maxLen }
					cnt[d]++
				}
				END{ for( i=1;i<=maxLen;i=i+1 ) printf "%d\t%d\n", i, cnt[i] }' \
			> $out
	else
		bamToBed -bedpe -i $src 2>&1 \
			| grep ^chr \
			| gawk 'BEGIN{ printf "fragLen\tCnt\n"; maxLen='$maxLen' }
				{ 
					if( $1=="." || $3=="." ) next
					# Distance between 5-prime ends
					if( $9=="+" ){
						d=$6-$2
					}else{
						d=$3-$5
					}
					if( d>maxLen ){ d=maxLen }
					cnt[d]++
				}
				END{ for( i=1;i<=maxLen;i=i+1 ) printf "%d\t%d\n", i, cnt[i] }' \
			> $out
	fi
elif [ $ext = "gz" ];then
	zcat $src \
		| grep ^chr \
		| gawk 'BEGIN{ printf "fragLen\tCnt\n"; maxLen='$maxLen' }
			{
				d=$3-$2
				if( d>maxLen ){ d=maxLen }
				cnt[d]++
			}
			END{ for( i=1;i<=maxLen;i=i+1 ) printf "%d\t%d\n", i, cnt[i] }' \
		> $out
else
	echo -e "Error: Invalid input file with extension \"$ext\"" >&@
	exit 1
fi
