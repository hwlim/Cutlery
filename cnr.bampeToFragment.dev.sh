#!/usr/bin/env bash


##############################################################################
# Written by Hee-Wooong Lim
# - Convert paired-end BAM file into fragment bed file by connecting two reads
# - Assume that reads are sorted by query name
#
	
source $MYBASHLIB/commonBash.sh
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [bam]
Description: Convert paired-end BAM file into fragment bed file by connecting two reads
	Assume that reads reads are sorted by query name
Options:
	-o <outFile>: output file, must be set.  default=<bamFile without path & extension>.frag.bed.gz
		e.g. ../input.bam -> input.frag.bed.gz
	-l <fragLen>: resize the fragment around the center. 0 for no resize. default=0
	-s: Sort by sort -k1,1 -k2,2n -k3,3n. default=Off
	-m: Memory size for sort, e.g. 10G. default=5G" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
des=NULL
finalLen=0
sortBed=FALSE
sortMem=5G
while getopts ":o:l:m:s" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		l)
			fragLen=$OPTARG
			;;
		m)
			sortMem=$OPTARG
			;;
		s)
			sortBam=TRUE
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

if [ "$des" != "NULL" ];then
	srcFile=`basename $src`
	prefix=${srcFile%.gz}
	prefix=${prefix%.bed}
	des=${prefix}.frag.bed.gz
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
