#!/usr/bin/env bash


###########################################3
# Written by Hee-Wooong Lim
# 
# Classify nucleosome-free reads and nucleosomal reads based on fragment length
#	1. Nucleosome-free reads  (NFR)
#	2. Mono-nucleosome reads (MON)
#	3. Nucleosomal reads (NUC): mono, di, tri nucleosomes
# For each category above, pair-connected files and pair-separate files are created
# => total 6 files are created
#

	
source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam]
Description: Split BAM file records into tree bed files
Options:
	-o <outPrefix>: Output prefix. required
	-f <samFlag>: SAM flag to include, default=0x2 (Properly paired)
	-F <samFlag>: SAM flag to exclude, default=0x400 (Duplicate)
	-l <finalLen>: final fixed size around the fragment center for *ctr.bed files, default=100
	-c <chromosome regex>: Regular expression for chromosome selection, default=^chr[0-9XY]+$
		For multiple patterns use regular expression, such as \"^chr[0-9XY]+$|^chrM$\" 
		NULL if not applicable or no filtering
Output:
	** Output files are not sorted **
	Read from nucleosome free regions
	- <outPrefix>.nfr.sep.bed: seprate read 1/2
	- <outPrefix>.nfr.con.bed: joined read 1/2 into single fragment
	- <outPrefix>.nfr.ctr.bed: joing read 1/2 into single fragment then resized to fixed length
	Read from nucleosomal regions
	- <outPrefix>.nuc.sep.bed: seprate read 1/2
	- <outPrefix>.nuc.con.bed: joined read 1/2 into single fragment
	- <outPrefix>.nuc.ctr.bed: joing read 1/2 into single fragment then resized to fixed length" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
outPrefix=NULL
flagInc=0x2
flagExc=0x400
finalLen=100
chrRegex='^chr[0-9XY]+$'
while getopts ":o:f:F:l:c:" opt; do
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

optStr=""
if [ "$flagInc" != "NULL" ];then
	optStr="-f $flagInc"
fi
if [ "$flagExc" != "NULL" ];then
	optStr="$optStr -F $flagExc"
fi

if [ "$chrRegex" == "NULL" ];then
	chrRegex="."
fi

printBAM(){
	local src=$1
	assertFileExist $src
	
	if [ "$optStr" == "" ];then
		bamToBed -bedpe -i $src
	else
		samtools view -u -b $optStr $src | bamToBed -bedpe
	fi
}

desNfr=${outPrefix}.nfr.sep.bed.gz
desNfrCon=${outPrefix}.nfr.con.bed.gz
desNfrCtr=${outPrefix}.nfr.ctr.bed.gz
desNuc=${outPrefix}.nuc.sep.bed.gz
desNucCon=${outPrefix}.nuc.con.bed.gz
desNucCtr=${outPrefix}.nuc.ctr.bed.gz

#######################################
## Start processing
echo -e "Splitting $src" >&2
echo -e "- SAM flag option = $optStr" >&2
echo -e "- Center fragment size = $finalLen" >&2
echo -e "- chrRegex = $chrRegex" >&2
echo -e "- Nucleosome free reads" >&2
echo -e "\t=> $desNfr" >&2
echo -e "\t=> $desNfrCon" >&2
echo -e "\t=> $desNfrCtr" >&2
echo -e "- Nucleosomal reads" >&2
echo -e "\t=> $desNuc" >&2
echo -e "\t=> $desNucCon" >&2
echo -e "\t=> $desNucCtr" >&2
#echo -e "TMPDIR = $TMPDIR" >&2

tmpNfr=${TMPDIR}/__temp__.$$.nfr
tmpNfrCon=${TMPDIR}/__temp__.$$.nfrCon
tmpNfrCtr=${TMPDIR}/__temp__.$$.nfrCtr
tmpNuc=${TMPDIR}/__temp__.$$.nuc
tmpNucCon=${TMPDIR}/__temp__.$$.nucCon
tmpNucCtr=${TMPDIR}/__temp__.$$.nucCtr

printBAM $src \
	| gawk 'BEGIN{
			tmpNfr="'$tmpNfr'"
			tmpNfrCon="'$tmpNfrCon'"
			tmpNfrCtr="'$tmpNfrCtr'"
			tmpNuc="'$tmpNuc'"
			tmpNucCon="'$tmpNucCon'"
			tmpNucCtr="'$tmpNucCtr'"

			printf "" > tmpNfr
			printf "" > tmpNfrCon
			printf "" > tmpNfrCtr
			printf "" > tmpNuc
			printf "" > tmpNucCon
			printf "" > tmpNucCtr
			hWid='$finalLen'/2
		}
		$1 ~ /'$chrRegex'/ {
			d=$6-$2
			c=($6+$2)/2
			if( d < 120 ){
				printf "%s\t%d\t%d\tNFR.%d_1\t0\t%s\n", $1,$2,$3,NR,$9 >> tmpNfr
				printf "%s\t%d\t%d\tNFR.%d_2\t0\t%s\n", $4,$5,$6,NR,$10 >> tmpNfr
				printf "%s\t%d\t%d\tNFR.%d\t0\t+\n", $1,$2,$6,NR >> tmpNfrCon
				printf "%s\t%d\t%d\tNFR.%d\t0\t+\n", $1,c-hWid,c+hWid,NR >> tmpNfrCtr
			}else if( d > 150 ){
				printf "%s\t%d\t%d\tNUC.%d_1\t0\t%s\n", $1,$2,$3,NR,$9 >> tmpNuc
				printf "%s\t%d\t%d\tNUC.%d_2\t0\t%s\n", $4,$5,$6,NR,$10 >> tmpNuc
				printf "%s\t%d\t%d\tNUC.%d\t0\t+\n", $1,$2,$6,NR >> tmpNucCon
				printf "%s\t%d\t%d\tNUC.%d\t0\t+\n", $1,c-hWid,c+hWid,NR >> tmpNucCtr
			}
		}'

gzip $tmpNfr
gzip $tmpNfrCon
gzip $tmpNfrCtr
gzip $tmpNuc
gzip $tmpNucCon
gzip $tmpNucCtr
mv ${tmpNfr}.gz $desNfr
mv ${tmpNfrCon}.gz $desNfrCon
mv ${tmpNfrCtr}.gz $desNfrCtr
mv ${tmpNuc}.gz $desNuc
mv ${tmpNucCon}.gz $desNucCon
mv ${tmpNucCtr}.gz $desNucCtr
