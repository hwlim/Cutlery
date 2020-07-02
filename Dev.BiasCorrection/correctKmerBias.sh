#!/usr/bin/env bash

###########################################################
# k-mer bias correction for bigwig file
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT
source $COMMON_LIB_BASE/commonBash.sh


function printUsage {
	echo -e "Usage: `basename $0` [options] <bigwig prefix>
Description:
	Perform k-mer bias correction for a given pair of bigWig files
Input:
	bigWig prefix assuming following two files
	- <prefix>.plus.bw
	- <prefix>.minus.bw
Options:
	-o <outPrefix>: outPrefix including path, required
	-l <k_up>: Upstream k-mer length, default=0
	-r <k_dn>: Downstream k-mer lengt includg 5'-end, default=0
		-----|--k_up--|--k_dn--|----
			       ================>
			       ^
			       5'-end
	-k : k-mer scaling file containing at least two columns in 1st/2nd columns without header, required
		k-mer / scaling factor (to multiply)
	-g <genome>: genome fasta file, required
	-s <chromSize>: chromosome size file, required
	-v : Verbose mode
Output:
	Bias-corrected bigWig files
	- <outPrefix>.plus.bw
	- <outPrefix>.minus.bw" >&2
}



###################################
## option and input file handling
k_up=0
k_dn=0
src_scaleFactor=NULL
genome=NULL
chromSize=NULL
verbose=FALSE
while getopts ":o:l:r:k:g:s:v" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
			;;
		l)
			k_up=$OPTARG
			;;
		r)
			k_dn=$OPTARG
			;;
		k)
			src_scaleFactor=$OPTARG
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

## Option & input validataion
bwPrefix=$1
bwPlus=${bwPrefix}.plus.bw
bwMinus=${bwPrefix}.minus.bw
assertFileExist $bwPlus $bwMinus

if [ "$genome" = "NULL" ];then
	echo -e "Error: genome fasta file (-g) must be specified" >&2
	exit 1
fi

assertFileExist $genome
if [ "$chromSize" = "NULL" ];then
	echo -e "Error: chrom.size file (-s) must be specified" >&2
	exit 1
fi
assertFileExist $chromSize

if [ "$src_scaleFactor" = "NULL" ];then
	echo -e "Error: scale factor file (-k) must be specified" >&2
	exit 1
fi
assertFileExist $src_scaleFactor

## k-mer length
if [ $k_up -eq 0 ] && [ $k_dn -eq 0 ];then
	echo -e "Both k_up and k_dn cannot be zero. Total k-mer length must be > 0" >&2
	exit 1
fi

## output validation
if [ "$outPrefix" = "NULL" ];then
	echo -e "Error: outPrefix (-o) must be specified" >&2
	exit 1
fi
desPlus=${outPrefix}.plus.bw
desMinus=${outPrefix}.minus.bw
if [ -f $desPlus ];then
	echo -e "Error: $desPlus already exists" >&2
	exit 1
fi
if [ -f $desMinus ];then
	echo -e "Error: $desMinus already exists" >&2
	exit 1
fi


## Start processing
if [ "$verbose" = "TRUE" ];then
	echo -e "==============================================" >&2
	echo -e "Bias correction" >&2
	echo -e "  bwPrefix = $bwPrefix" >&2
	echo -e "  K-mer = $k_up / $k_dn" >&2
	echo -e "  genome = $genome" >&2
	echo -e "  desPrefix = $desPrefix" >&2
	echo -e "  scaleFactor = $scaleFactor" >&2
	echo -e "==============================================" >&2
fi

desDir=`dirname ${outPrefix}`
mkdir -p $desDir



correctBias(){
	local src=$1
	local des=$2
	local strand=$3

	## Temporary files
	tmpKmerBed=${TMPDIR}/__temp__.$$.kmerBed.txt
	tmpCorrected=${TMPDIR}/__temp__.$$.corrected.txt
	tmpBg=${TMPDIR}/__temp__.$$.bg
	tmpBw=${TMPDIR}/__temp__.$$.bw


	## k-mer region extraction
	## 6-column bed format output
	## 5th column is bigWig signal value
	## bedtools cannot handle chr:0-0, so explicitly filtering out
	[ "$verbose" = "TRUE" ] && echo -e "  - Making $bed_kmer" >&2
	if [ "$strand" = "plus" ];then
		gawkStr='{ if($2<1) next; for( i=$2; i<$3; i=i+1 ) printf "%s\t%d\t%d\t%s\t0\t+\n", $1, i, i, $4 }'
	elif [ "$strand" = "minus" ];then
		gawkStr='{ if($2<1) next; for( i=$2+1; i<=$3; i=i+1 ) printf "%s\t%d\t%d\t%s\t0\t-\n", $1, i, i, $4 }'
	else
		echo -e "Error: Invalid strand: $strand" >&2
	fi

	bigWigToBedGraph $src /dev/stdout \
		| gawk "$gawkStr" \
		| slopBed -g $chromSize -s -l ${k_up} -r ${k_dn} \
		| gawk '$3-$2=='$kmer'' \
		> $tmpKmerBed
	#if [ "$strand" = "plus" ];then
	#	bigWigToBedGraph $src /dev/stdout \
	#		| gawk '{ for( i=$2; i<$3; i=i+1 ) printf "%s\t%d\t%d\t%s\t0\t+\n", $1, i, i, $4 }' \
	#		| gawk '$2!=0' \
	#		| slopBed -g $chromSize -s -l ${k_up} -r ${k_dn} \
	#		| gawk '$3-$2=='$kmer'' \
	#		> $tmp
	#elif [ "$strand" = "minus" ];then
	#	bigWigToBedGraph $src /dev/stdout \
	#		| gawk '{ for( i=$2+1; i<=$3; i=i+1 ) printf "%s\t%d\t%d\t%s\t0\t-\n", $1, i, i, $4 }'\
	#		| gawk '$2!=0' \
	#		| slopBed -g $chromSize -s -l ${k_up} -r ${k_dn} \
	#		| gawk '$3-$2=='$kmer'' \
	#		> $tmp
	#else
	#	echo -e "Error: Invalid strand: $strand" >&2
	#fi


	## k-mer extraction & signal correction
	[ "$verbose" = "TRUE" ] && echo -e "  - Correcting values" >&2
	bedtools getfasta -bed $tmpKmerBed -fi $genome -fo /dev/stdout -name -tab -s \
		| sed 's/::/\t/' \
		| cut -f 1,3 \
		| tr '[a-z]' '[A-Z]' \
		| gawk 'BEGIN{  while(getline < "'${src_scaleFactor}'"){ scaleDic[$1]=$2 }  }
			{
				printf "%s\t%.3f\n", $2, $1*scaleDic[$2]
			}' \
		> $tmpCorrected
	

	## Reconstraction of a corrected bedGraph file
	[ "$verbose" = "TRUE" ] && echo -e "  Creating bedGraph file" >&2
	if [ "$strand" = "plus" ];then
		gawkStr='{ c=$2+'${k_up}'; printf "%s\t%d\t%d\t%s\n", $1,c,c+1,$8 }'
	elif [ "$strand" = "minus" ];then
		gawkStr='{ c=$3-'${k_up}'; printf "%s\t%d\t%d\t%s\n", $1,c-1,c,$8 }'
	else
		echo -e "Error: Invalid strand: $strand" >&2
	fi
	paste ${tmpKmerBed} ${kmerCorrected} \
		| gawk "$gawkStr" \
		> $tmpBg
	rm $tmpKmerBed
	rm $kmerCorrected

	#echo -e "  - Making $bg_corrected" >&2
	#if [ "$strand" = "plus" ];then
	#	paste ${tmpKmerBed} ${kmerCorrected} \
	#		| gawk '{ c=$2+'${k_up}'; printf "%s\t%d\t%d\t%s\n", $1,c,c+1,$8 }' \
	#		> $tmp
	#elif [ "$strand" = "minus" ];then
	#	paste ${tmpKmerBed} ${kmerCorrected} \
	#		| gawk '{ c=$3-'${k_up}'; printf "%s\t%d\t%d\t%s\n", $1,c-1,c,$8 }' \
	#		> $tmp
	#else
	#	echo -e "Error: Invalid strand: $strand" >&2
	#fi
	#mv $tmp $bg_corrected

	## bedGraph to bigWig
	[ "$verbose" = "TRUE" ] && echo -e "  - Making $bw_corrected" >&2
	bedGraphToBigWig $tmpBg ${chromSize} $tmpBw
	mv $tmpBw $des
	rm $tmpBg
}


correctBias $bwPlus ${desPlus} plus
correctBias $bwMinus ${desMinus} minus

