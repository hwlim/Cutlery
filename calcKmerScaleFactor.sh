#!/usr/bin/env bash

###########################################################
# Scale factor calculation using k-mer count from ChIP-exo data & genomewide count
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT
source $COMMON_LIB_BASE/commonBash.sh


function printUsage {
	echo -e "Usage: `basename $0` <k-mer count file>
Description: Calculate scaling factor based on given kmer counts
Options:
	-g: genomic k-mer counts file, two columns (kmer / count)
	-p: pseudo count, default=1
	-v: Verbose mode
Input:
	Two column text file containing:
	1. k-mer sequence (in all upper case)
	2. k-mer count
	without header
Output:
	Four column output to STDOUT without header
	1. k-mer
	2. scale factor (for multiply)
	3. sample k-mer count
	4. genomic k-mer count" >&2
}



###################################
## option and input file handling
src_kmer_genome=NULL
pseudo=1
verbose=FALSE
while getopts ":g:p:v" opt; do
	case $opt in
		g)
			src_kmer_genome=$OPTARG
			;;
		p)
			pseudo=$OPTARG
			;;
		v)
			verbose=TRUE
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
if [ $# -lt 1 ];then
	printUsage
	exit 1
fi

src_kmer=$1
assertFileExist $src_kmer $src_kmer_genome

## k-mer length validation
set +o pipefail
kmer_len_sample=`head -n 1 $src_kmer | gawk '{ printf "%d", length($1) }'`
kmer_len_genome=`head -n 1 $src_kmer_genome | gawk '{ printf "%d", length($1) }'`
set -o pipefail
if [ $kmer_len_sample -ne $kmer_len_genome ];then
	echo -e "Error: k-mer length does not match between input file vs genome" >&2
	exit 1
fi

###### Input data format
## src_kmer file contents	: kmer / plus count / minus count / total count
## src_kmer_genome file contents: kmer / count

## Total k-mer counts in genome
ttcGenome=`cat $src_kmer_genome | gawk 'BEGIN{s=0}{ s = s + $2 }END{ printf "%d", s }'`

## Number of k-mers
N_kmer=`cat $src_kmer_genome | grep -v ^$ | wc -l`

## Total k-mer counts from sample; pseudo counts are added in calculating total read counts except for k-mers containing 'N'
ttc=`cat $src_kmer | grep -v -e N | gawk 'BEGIN{s=0}{ s=s+$4 }END{ printf "%d", s + '${N_kmer}'*'${pseudo}' }'`

if [ "${verbose}" = "TRUE" ];then
	echo -e "================================================">&2
	echo -e "Calculating k-mer scale factors" >&2
	echo -e "  - Total read counts = $ttc" >&2
	echo -e "  - k-mer count = $src_kmer" >&2
	echo -e "  - Genomic k-mer count = $src_kmer_genome" >&2
	echo -e "  - k-mer length = $kmer_len_sample" >&2
	echo -e "  - pseudo count = $pseudo" >&2
fi

cat $src_kmer_genome \
	| gawk 'BEGIN{
			while(getline < "'${src_kmer}'"){ cntDic[$1]=$4 }
			ttc='$ttc'; ttcGenome='$ttcGenome'; pseudo='$pseudo';
		}
		{
			kmer=$1
			cntGenome = $2
			if( kmer in cntDic ){
				cnt = cntDic[kmer] + pseudo
			}else{
				cnt = pseudo
			}

			ratioSample = cnt / ttc
			ratioGenome = cntGenome / ttcGenome
			scaleFactor = ratioGenome / ratioSample
			printf "%s\t%.6f\t%d\t%d\n", kmer, scaleFactor, cnt, cntGenome
		}'

