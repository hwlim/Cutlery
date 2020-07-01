#!/usr/bin/env bash

###########################################################
# k-mer bias correction for bigwig file
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT
source $COMMON_LIB_BASE/commonBash.sh


function printUsage {
	echo -e "Usage: `basename $0` [options] <bigwig prefix>
Description:
	Perform k-mer bias correction for a given pair of bigWig files
	-----|--k_up--|--k_dn--|----
	               ================>
	               ^
	               5'-end
Input:
	bigWig prefix assuming following two files
	- <prefix>.plus.bw
	- <prefix>.minus.bw
Options:
	-l <k_up>: Upstream k-mer length, default=0
	-r <k_dn>: Downstream k-mer lengt includg 5'-end, default=0
	-k : k-mer scaling file containing at least two columns in 1st/2nd columns
		k-mer / scaling factor (to multiply)
	-g <genome>: genome fasta file, required
	-s <chromSize>: chromosome size file, required
	-v : Verbose mode
Output:
	Two column output (k-mer and count) in STDOUT
Output:
	Four column output to STDOUT without header
	1. k-mer
	2. scale factor (for multiply)
	3. sample k-mer count
	4. genomic k-mer count" >&2
}



###################################
## option and input file handling
k_up=0
k_dn=0
src_scaling=NULL
genome=NULL
chromSize=NULL
verbose=FALSE
while getopts ":l:r:k:g:s:v" opt; do
	case $opt in
		l)
			k_up=$OPTARG
			;;
		r)
			k_dn=$OPTARG
			;;
		k)
			src_scaling=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		s)
			chromSize=$OPTARG
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

bwPrefix=$1
bwPlus=${bwPrefix}.plus.bw
bwMinus=${bwPrefix}.minus.bw

assertFileExist $genome $chromSize
assertFileExist $bwPlus $bwMinus

echo -e "==============================================" >&2
echo -e "Bias correction" >&2
echo -e "  bwPrefix = $bwPrefix" >&2
echo -e "  K-mer = $k_up / $k_dn" >&2
echo -e "  genome = $genome" >&2
echo -e "  desPrefix = $desPrefix" >&2
echo -e "  scaleFactor = $scaleFactor" >&2
echo -e "==============================================" >&2

tmp=__temp__.$$.txt

if [ "$bwPrefix" = "NULL" ];then
	echo -e "bwPrefix is NULL, skip" >&2
	exit 0
fi

correctBiasExo(){
	src=$1
	prefix=$2
	strand=$3

	echo -e ">> Processing $src" >&2

	bed_kmer=${prefix}.${strand}.bed
	kmer_corrected=${prefix}.${strand}.txt
	bg_corrected=${prefix}.${strand}.bg
	bw_corrected=${prefix}.${strand}.bw

	if [ -f $bw_corrected ];then
		echo -e "  - $bw_corrected already exist, pass" >&2
		return
	fi

	## k-mer region extraction
	## 6-column bed format output
	## 5th column is bigWig signal value
	if [ -f $bed_kmer ] || [ -f ${bed_kmer}.gz ];then
		echo -e "  - $bed_kmer already exist, pass" >&2
	else
		## bedtools cannot handle chr:0-0, so explicitly filtering out
		echo -e "  - Making $bed_kmer" >&2
		if [ "$strand" = "plus" ];then
			bigWigToBedGraph $src /dev/stdout \
				| gawk '{ for( i=$2; i<$3; i=i+1 ) printf "%s\t%d\t%d\t%s\t0\t+\n", $1, i, i, $4 }' \
				| gawk '$2!=0' \
				| slopBed -g $chromSize -s -l ${k_up} -r ${k_dn} \
				| gawk '$3-$2=='$kmer'' \
				> $tmp
		elif [ "$strand" = "minus" ];then
			bigWigToBedGraph $src /dev/stdout \
				| gawk '{ for( i=$2+1; i<=$3; i=i+1 ) printf "%s\t%d\t%d\t%s\t0\t-\n", $1, i, i, $4 }'\
				| gawk '$2!=0' \
				| slopBed -g $chromSize -s -l ${k_up} -r ${k_dn} \
				| gawk '$3-$2=='$kmer'' \
				> $tmp
		else
			echo -e "Error: Invalid strand: $strand" >&2
		fi
		mv $tmp ${bed_kmer}
	fi


	## k-mer extraction & signal correction
	if [ -f $kmer_corrected ];then
		echo -e "  - $kmer_corrected already exist, pass" >&2
	else
		echo -e "  - Making $kmer_corrected" >&2
		bedtools getfasta -bed $bed_kmer -fi $genomeFa -fo /dev/stdout -name -tab -s \
			| sed 's/::/\t/' \
			| cut -f 1,3 \
			| tr '[a-z]' '[A-Z]' \
			| gawk 'BEGIN{  while(getline < "'${scaleFactor}'"){ scaleDic[$1]=$2 }  }
				{
					printf "%s\t%.3f\n", $2, $1*scaleDic[$2]
				}' \
			> $tmp
		mv $tmp $kmer_corrected
	fi

	## Reconstraction of a corrected bedGraph file
	if [ -f $bg_corrected ] || [ -f ${bg_corrected}.gz ];then
		echo -e "  - $bg_corrected already exist, pass" >&2
	else
		echo -e "  - Making $bg_corrected" >&2
		if [ "$strand" = "plus" ];then
			paste ${bed_kmer} ${kmer_corrected} \
				| gawk '{ c=$2+'${k_up}'; printf "%s\t%d\t%d\t%s\n", $1,c,c+1,$8 }' \
				> $tmp
		elif [ "$strand" = "minus" ];then
			paste ${bed_kmer} ${kmer_corrected} \
				| gawk '{ c=$3-'${k_up}'; printf "%s\t%d\t%d\t%s\n", $1,c-1,c,$8 }' \
				> $tmp
		else
			echo -e "Error: Invalid strand: $strand" >&2
		fi
		mv $tmp $bg_corrected
	fi

	## bedGraph to bigWig
	if [ -f $bw_corrected ];then
		echo -e "  - $bw_corrected already exist, pass" >&2
	else
		echo -e "  - Making $bw_corrected" >&2
		if [ -f ${bg_corrected}.gz ];then
			gunzip -v ${bg_corrected}.gz
		fi
		bedGraphToBigWig $bg_corrected ${chromSize} $tmp
		mv $tmp $bw_corrected
		gzip -v $bg_corrected
	fi
}

bwPlus=${bwPrefix}.plus.bw
bwMinus=${bwPrefix}.minus.bw
assertFileExist $bwPlus $bwMinus

correctBiasExo $bwPlus ${desPrefix} plus
correctBiasExo $bwMinus ${desPrefix} minus

ln -f -s ${bwPlus} ${desPrefixRaw}.plus.bw
ln -f -s ${bwMinus} ${desPrefixRaw}.minus.bw
