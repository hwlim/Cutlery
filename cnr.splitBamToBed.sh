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

#########################
# Not to considering for future revision
# - Currently, read1/read2 in a bed file are not actual Read1/Read2 from the origianl bam file
#	Instead, read1 is left side read and read2 is right side read after alignment.
# 	We may need a minor revision to preserve Read1/Read2 information in the final bed file
# - Protrusion handling
#	Currently, protrusion is appropriately handled in making fragment bed file
#	However, it is not being handled in creating seprate (sep) bed file
#	Each read (R1/R2) files start/end coordinate should be correced in consistent with fragments



	
source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bam]
Description: Split BAM file records into three bed files
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
	All reads
	- <outPrefix>.all.sep.bed: seprate read 1/2
	- <outPrefix>.all.con.bed: joined read 1/2 into single fragment
	- <outPrefix>.all.ctr.bed: joing read 1/2 into single fragment then resized to fixed length
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

desAll=${outPrefix}.all.sep.bed.gz
desAllCon=${outPrefix}.all.con.bed.gz
desAllCtr=${outPrefix}.all.ctr.bed.gz

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
echo -e "- All reads" >&2
echo -e "\t=> $desAll" >&2
echo -e "\t=> $desAllCon" >&2
echo -e "\t=> $desAllCtr" >&2
echo -e "- Nucleosome free reads" >&2
echo -e "\t=> $desNfr" >&2
echo -e "\t=> $desNfrCon" >&2
echo -e "\t=> $desNfrCtr" >&2
echo -e "- Nucleosomal reads" >&2
echo -e "\t=> $desNuc" >&2
echo -e "\t=> $desNucCon" >&2
echo -e "\t=> $desNucCtr" >&2
#echo -e "TMPDIR = $TMPDIR" >&2

tmpAll=${TMPDIR}/__temp__.$$.all
tmpAllCon=${TMPDIR}/__temp__.$$.allCon
tmpAllCtr=${TMPDIR}/__temp__.$$.allCtr
tmpNfr=${TMPDIR}/__temp__.$$.nfr
tmpNfrCon=${TMPDIR}/__temp__.$$.nfrCon
tmpNfrCtr=${TMPDIR}/__temp__.$$.nfrCtr
tmpNuc=${TMPDIR}/__temp__.$$.nuc
tmpNucCon=${TMPDIR}/__temp__.$$.nucCon
tmpNucCtr=${TMPDIR}/__temp__.$$.nucCtr

printBAM $src \
	| gawk 'BEGIN{
			tmpAll="'$tmpAll'"
			tmpAllCon="'$tmpAllCon'"
			tmpAllCtr="'$tmpAllCtr'"
			tmpNfr="'$tmpNfr'"
			tmpNfrCon="'$tmpNfrCon'"
			tmpNfrCtr="'$tmpNfrCtr'"
			tmpNuc="'$tmpNuc'"
			tmpNucCon="'$tmpNucCon'"
			tmpNucCtr="'$tmpNucCtr'"

			printf "" > tmpAll
			printf "" > tmpAllCon
			printf "" > tmpAllCtr
			printf "" > tmpNfr
			printf "" > tmpNfrCon
			printf "" > tmpNfrCtr
			printf "" > tmpNuc
			printf "" > tmpNucCon
			printf "" > tmpNucCtr
			hWid='$finalLen'/2
		}
		$1 ~ /'$chrRegex'/ {
			if( $9=="+" ){
				d=$6-$2
				s=$2; e=$6
			}else{
				d=$3-$5
				s=$5; e=$3
			}
			c=(s+e)/2			
			#d=$6-$2
			#c=($6+$2)/2
			printf "%s\t%d\t%d\tAll.%d/1\t0\t%s\n", $1,$2,$3,NR,$9 >> tmpAll
			printf "%s\t%d\t%d\tAll.%d/2\t0\t%s\n", $4,$5,$6,NR,$10 >> tmpAll
			printf "%s\t%d\t%d\tAll.%d\t0\t+\n", $1,s,e,NR >> tmpAllCon
			#printf "%s\t%d\t%d\tAll.%d\t0\t+\n", $1,$2,$6,NR >> tmpAllCon
			printf "%s\t%d\t%d\tAll.%d\t0\t+\n", $1,c-hWid,c+hWid,NR >> tmpAllCtr
			if( d < 120 ){
				printf "%s\t%d\t%d\tNFR.%d/1\t0\t%s\n", $1,$2,$3,NR,$9 >> tmpNfr
				printf "%s\t%d\t%d\tNFR.%d/2\t0\t%s\n", $4,$5,$6,NR,$10 >> tmpNfr
				printf "%s\t%d\t%d\tNFR.%d\t0\t+\n", $1,s,e,NR >> tmpNfrCon
				#printf "%s\t%d\t%d\tNFR.%d\t0\t+\n", $1,$2,$6,NR >> tmpNfrCon
				printf "%s\t%d\t%d\tNFR.%d\t0\t+\n", $1,c-hWid,c+hWid,NR >> tmpNfrCtr
			}else if( d > 150 ){
				printf "%s\t%d\t%d\tNUC.%d/1\t0\t%s\n", $1,$2,$3,NR,$9 >> tmpNuc
				printf "%s\t%d\t%d\tNUC.%d/2\t0\t%s\n", $4,$5,$6,NR,$10 >> tmpNuc
				printf "%s\t%d\t%d\tNUC.%d\t0\t+\n", $1,s,e,NR >> tmpNucCon
				#printf "%s\t%d\t%d\tNUC.%d\t0\t+\n", $1,$2,$6,NR >> tmpNucCon
				printf "%s\t%d\t%d\tNUC.%d\t0\t+\n", $1,c-hWid,c+hWid,NR >> tmpNucCtr
			}
		}'

gzip $tmpAll
gzip $tmpAllCon
gzip $tmpAllCtr
gzip $tmpNfr
gzip $tmpNfrCon
gzip $tmpNfrCtr
gzip $tmpNuc
gzip $tmpNucCon
gzip $tmpNucCtr
mv ${tmpAll}.gz $desAll
mv ${tmpAllCon}.gz $desAllCon
mv ${tmpAllCtr}.gz $desAllCtr
mv ${tmpNfr}.gz $desNfr
mv ${tmpNfrCon}.gz $desNfrCon
mv ${tmpNfrCtr}.gz $desNfrCtr
mv ${tmpNuc}.gz $desNuc
mv ${tmpNucCon}.gz $desNucCon
mv ${tmpNucCtr}.gz $desNucCtr
