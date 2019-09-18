#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam or bed(.gz)]
Description:
	Create stranded bigWig files using input alignment bam or bed(.gz) file
Input:
	bam or bed(.gz) file
Options:
	-o <outPrefix>: Output file prefix including path, default=<src file name>
		<outPrefix>.plus.bw
		<outPrefix>.minus.bw
	-g <chrom.size>: chromosome size file, required
	-l <fragLen>: fragment length. Each alignment is resized to this length toward 3'-end. default=0 (as it is)
	-m <sortMem>: memory for sorting, default=5G
	-c <chromosome regex>: Regular expression for chromosome selection, default=^chr[0-9XY]+$
		For multiple patterns use regular expression, such as \"^chr[0-9XY]+$|^chrM$\" 
		NULL if not applicable or no filtering" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desPrefix=""
chrom=""
fragLen=0
sortMem=5G
chrRegex='^chr[0-9XY]+$'
while getopts ":o:g:l:m:c:" opt; do
	case $opt in
		o)
			desPrefix=$OPTARG
			;;
		g)
			chrom=$OPTARG
			;;
		l)
			fragLen=$OPTARG
			;;
		m)
			sortMem=$OPTARG
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


if [ "$chrom" == "" ];then
	echo -e "Error: chrom.size file (-g) must be specified" >&2
	exit 1
fi

assertFileExist $src $chrom


if [ "$desPrefix" == "" ];then
	srcFile=`basename $src`
	desPrefix=`echo $srcFile | sed 's/.bam$\|.bed.*$//'`
fi
desDir=`dirname $desPrefix`
mkdir -p $desDir


echo -e "Making stranded bigWig files" >&2
echo -e "  src = $src" >&2
echo -e "  chrom = $chrom" >&2
echo -e "  fragLen = $fragLen" >&2
echo -e "  sortMem = $sortMem" >&2
echo -e "  chrRegex = $chrRegex" >&2

printAlign(){
	if [ -d "$src" ];then
		## Homer TagDir
		echo -e "$src is a directory; assume it as a homer tag directory" >&2
		tagDir2bed.pl -separate $src
	elif [ -f "$src" ];then
		## BED(.gz) or BAM file
		ext=${src##*.}
		if [ "$ext" == "gz" ];then
			zcat $src
		elif [ "$ext" == "bed" ];then
			cat $src
		elif [ "$ext" == "bam" ];then
			bamToBed -i $src
		else
			echo -e "Error: $src is unknown format" >&2
			exit 1
		fi
	else
		echo -e "Error: $src does not exist" >&2
		exit 1
	fi
}


tmpChrom=${TMPDIR}/$$.chrom.size
bed=${TMPDIR}/$$.bed.gz
tmpBg_plus=${TMPDIR}/$$.plus.bedGraph
tmpBg_minus=${TMPDIR}/$$.minus.bedGraph
tmpBw_plus=${TMPDIR}/$$.plus.bw
tmpBw_minus=${TMPDIR}/$$.minus.bw
bw_plus=${name}.plus.bw
bw_minus=${name}.minus.bw

echo $chrom \
	| gawk '$1 ~ /'$chrRegex'/' \
	| sort -k1,1 \
	> $tmpChrom

## Resized alignment bed file
echo -e "1) Making (resized) bed file: fragLen=${fragLen}" >&2
if [ ${fragLen} -gt 0 ];then
	printAlign ${src} \
		| gawk '{
				if( $6 == "+" ){
					s = $2
					e = s + '$fragLen'
				}else{
					e = $3
					s = e - '$fragLen'
				}
				printf "%s\t%d\t%d\t.\t0\t%s\n", $1, s, e, $6
			}' \
		| sort -S $sortMem -k1,1 -k2,2n -k3,3n \
		> $bed
else
	printAlign ${src} \
		| sort -S $sortMem -k1,1 -k2,2n -k3,3n \
		> $bed
fi


## Total tag count
ttc=`zcat $bed | wc -l`

## BedGraph file
echo -e "2) Making bedGraph files\n\t TotalTags: $ttc" >&2
echo -e "\tStep 1: Generating Plus: $tmpBg_plus" >&2
genomeCoverageBed -bg -i ${bed} -strand + -g $tmpChrom \
	| gawk 'BEGIN{ttc='$ttc'}
		{
			printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4*1000000/ttc;
		}' \
	> ${tmpBg_plus}
echo -e "\tStep 2: Generating Minus: $tmpBg_minus" >&2
genomeCoverageBed -bg -i ${bed} -strand - -g $tmpChrom \
	| gawk 'BEGIN{ttc='$ttc'}
		{
			printf "%s\t%s\t%s\t-%.5f\n", $1,$2,$3,$4*1000000/ttc;
		}' \
	> ${tmpBg_minus}


## BigWig file
isFileExist $tmpBw_plus $tmpBw_minus
echo -e "3) Making bigWig files" >&2
echo -e "\tStep 1: Generating Plus: $bw_plus" >&2
bedGraphToBigWig ${tmpBg_plus} ${tmpChrom} ${tmpBw_plus}
echo -e "\tStep 2: Generating Minus: $bw_minus" >&2
bedGraphToBigWig ${tmpBg_minus} ${tmpChrom} ${tmpBw_minus}

mv ${tmpBw_plus} ${bw_plus}
mv ${tmpBw_minus} ${bw_minus}
