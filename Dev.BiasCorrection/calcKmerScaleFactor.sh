#!/usr/bin/env bash

###########################################################
# Scale factor calculation using k-mer count from ChIP-exo data & genomewide count
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT
source $MYBASHLIB/commonBash.sh
# kmerCnt
# kmerCntGenome
# scaleFactor



function printUsage {
	echo -e "Usage: `basename $0` <kmerCntGenome> <kmerCnt>" >&2
	echo -e "Description: Calculate scaling factor based on given kmer counts" >&2
	echo -e "Options:" >&2
	echo -e "\t-v : Verbose mode" >&2
	echo -e "Output:" >&2
	echo -e "\tFour column output (kmer / scale factor / exo-Kmer / genome-Kmer) in STDOUT" >&2
}



###################################
## option and input file handling
verbose=FALSE
while getopts ":v" opt; do
	case $opt in
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
if [ $# -lt 2 ];then
	printUsage
	exit 1
fi

kmerCntGenome=$1
kmerCnt=$2

assertFileExist $kmerCnt $kmerCntGenome


pseudo=1

ttcGen=`cat $kmerCntGenome | gawk 'BEGIN{s=0}{ s=s+$3 }END{ printf "%d", s }'`
N_kmer=`cat $kmerCntGenome | wc -l`
ttcExo=`cat $kmerCnt | grep -v N | gawk 'BEGIN{s=0}{ s=s+$2 }END{ printf "%d", s + '${N_kmer}'*'${pseudo}' }'`
# Pseudo-count 1 is added for each frequency value of k-mer

if [ "${verbose}" = "TRUE" ];then
	echo -e "==============================================">&2
	echo -e "Calculating scale factors by K-mer Total Counts" >&2
	echo -e "  Exo = $ttcExo" >&2
	echo -e "  Genome = $ttcGen" >&2
fi

## kmerCntGenome file contents
# 1	AAAAAA	5328209	5076870
cat $kmerCntGenome \
	| gawk 'BEGIN{
			while(getline < "'${kmerCnt}'"){ cntDic[$1]=$2 }
			ttcExo='$ttcExo'; ttcGen='$ttcGen';
		}
		{
			if( $2 in cntDic ){
				cnt = cntDic[$2] + '${pseudo}'
			}else{
				cnt = '${pseudo}'
			}
			exo=cnt/ttcExo
			gen=$3/ttcGen
			factor= gen / exo
			printf "%s\t%.6f\t%d\t%d\n", $2, factor, cnt, $3

		}'

