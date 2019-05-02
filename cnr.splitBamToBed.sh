#!/usr/bin/env bash


###########################################3
# ATAC-seq tools
# Written by Hee-Wooong Lim
# 
# Classify nucleosome-free reads and nucleosomal reads based on fragment length
#	1. Nucleosome-free reads  (NFR)
#	2. Mono-nucleosome reads (MON)
#	3. Nucleosomal reads (NUC): mono, di, tri nucleosomes
# For each category above, pair-connected files and pair-separate files are created
# => total 6 files are created
#
	
source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [inputFiles] ..." >&2
	echo -e "Description: Split BAM file records into tree bed files" >&2
	echo -e "Options:" >&2
        echo -e "\t-o <outDir>: Destination directory, default=<same with the src file>" >&2
	echo -e "\t-f <samFlag>: SAM flag to include, default=0x2 (Properly paired)" >&2
	echo -e "\t-F <samFlag>: SAM flag to exclude, default=0x400 (Duplicate)" >&2
	echo -e "\t-l <finalLen>: final fixed size around the fragment center for *ctr.bed files, default=100" >&2
	echo -e "Output:" >&2
	echo -e "\tRead from nucleosome free regions" >&2
	echo -e "\t- <filename>.nfr.sep.bed: seprate read 1 & 2" >&2
	echo -e "\t- <filename>.nfr.con.bed: joined read 1 & 2 into single fragment" >&2
	echo -e "\t- <filename>.nfr.ctr.bed: joing read 1 & 2 into single fragment then resized to fixed length" >&2
	echo -e "\tRead from nucleosomal regions" >&2
	echo -e "\t- <filename>.nuc.sep.bed: seprate read 1 & 2" >&2
	echo -e "\t- <filename>.nuc.con.bed: joined read 1 & 2 into single fragment" >&2
	echo -e "\t- <filename>.nuc.ctr.bed: joing read 1 & 2 into single fragment then resized to fixed length" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desBase=NULL
flagInc=0x2
flagExc=0x400
finalLen=100
while getopts ":o:f:F:l:" opt; do
	case $opt in
		o)
			desBase=$OPTARG
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

assertFileExist $@


###################################
## main code

if [ "$desBase" != "NULL" ];then
	mkdir -p $desBase
fi

optStr=""
if [ "$flagInc" != "NULL" ];then
	optStr="-f $flagInc"
fi
if [ "$flagExc" != "NULL" ];then
	optStr="$optStr -F $flagExc"
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

for src in $@
do
	srcFile=`basename $src`
	if [ "$desBase" = "NULL" ];then
		desDir=`dirname $src`
	else
		desDir=$desBase
	fi
	desNfr=${desDir}/${srcFile%.bam}.nfr.sep.bed
	desNfrCon=${desDir}/${srcFile%.bam}.nfr.con.bed
	desNfrCtr=${desDir}/${srcFile%.bam}.nfr.ctr.bed
	desNuc=${desDir}/${srcFile%.bam}.nuc.sep.bed
	desNucCon=${desDir}/${srcFile%.bam}.nuc.con.bed
	desNucCtr=${desDir}/${srcFile%.bam}.nuc.ctr.bed

	echo -e "Splitting $src" >&2
	echo -e "  SAM flag option = $optStr" >&2
	echo -e "  Nucleosome free reads" >&2
	echo -e "\t=> $desNfr" >&2
	echo -e "\t=> $desNfrCon" >&2
	echo -e "\t=> $desNfrCtr" >&2
	echo -e "  Nucleosomal reads" >&2
	echo -e "\t=> $desNuc" >&2
	echo -e "\t=> $desNucCon" >&2
	echo -e "\t=> $desNucCtr" >&2

	tmpNfr=__temp__.$$.nfr
	tmpNfrCon=__temp__.$$.nfrCon
	tmpNfrCtr=__temp__.$$.nfrCtr
	tmpNuc=__temp__.$$.nuc
	tmpNucCon=__temp__.$$.nucCon
	tmpNucCtr=__temp__.$$.nucCtr

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
			{
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

	mv $tmpNfr $desNfr
	mv $tmpNfrCon $desNfrCon
	mv $tmpNfrCtr $desNfrCtr
	mv $tmpNuc $desNuc
	mv $tmpNucCon $desNucCon
	mv $tmpNucCtr $desNucCtr
	echo -e "" >&2
done
