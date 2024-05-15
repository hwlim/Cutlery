#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 
# Split input fragment file into two portions by size
#	- NFR (< 120bp) and NUC (>150bp)


	
source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [fragment.bed.gz]
Description: Split BAM file records into three bed files
Options:
	-o <outPrefix>: Output prefix. required
	-l <finalLen>: final fixed size around the fragment center, default=100. (Set 0 to keep the original size)
	-c <chromosome regex>: Regular expression for chromosome selection, default=^chr[0-9XY]+$
		For multiple patterns use regular expression, such as \"^chr[0-9XY]+$|^chrM$\" 
		NULL if not applicable or no filtering
Output:
	** Note: output bed files are not sorted **
	- <outPrefix>.nfr.bed.gz: NFR fragments
	- <outPrefix>.nuc.bed.gz: NUC fragments" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
outPrefix=NULL
finalLen=100
chrRegex='^chr[0-9XY]+$'
csorted=FALSE
while getopts ":o:l:c:s" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
			;;
		f)
			flagInc=$OPTARG
			;;
		F)
			flagExc=$OPTARG
			;;
		l)
			finalLen=$OPTARG
			;;
		c)
			chrRegex=$OPTARG
			;;
		s)
			csorted=TRUE
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


###################################
## main code

if [ "$outPrefix" == "NULL" ];then
	echo -e "Error: outPrefix must be set" >&2
	exit 1
fi

desDir=`dirname $outPrefix`
mkdir -p $desDir

if [ "$chrRegex" == "NULL" ];then
	chrRegex="."
fi

desNfr=${outPrefix}.nfr.bed.gz
desNuc=${outPrefix}.nuc.bed.gz


#######################################
## Start processing
echo -e "Splitting $src" >&2
echo -e "- Center fragment size = $finalLen" >&2
echo -e "- chrRegex = $chrRegex" >&2
echo -e "- Output" >&2
echo -e "\t=> $desNfr" >&2
echo -e "\t=> $desNuc" >&2

tmpNfr=${TMPDIR}/__temp__.$$.nfr
tmpNuc=${TMPDIR}/__temp__.$$.nuc

zcat $src \
	| gawk 'BEGIN{
			tmpNfr="'$tmpNfr'"
			tmpNuc="'$tmpNuc'"

			printf "" > tmpNfr
			printf "" > tmpNuc
			hWid='$finalLen'/2
		}
		$1 ~ /'$chrRegex'/ {
			s = $2
			e = $3
			d = e - s
			c=(s+e)/2	

			if(hWid>0){
				s = c - hWid
				e = c + hWid
			}

			if( d < 120 ){
				printf "%s\t%d\t%d\tNFR.%d\t0\t+\n", $1,s,e,NR >> tmpNfr
			}else if( d > 150 ){
				printf "%s\t%d\t%d\tNUC.%d\t0\t+\n", $1,s,e,NR >> tmpNuc
			}
		}'

gzip $tmpNfr
gzip $tmpNuc
mv ${tmpNfr}.gz $desNfr
mv ${tmpNuc}.gz $desNuc
