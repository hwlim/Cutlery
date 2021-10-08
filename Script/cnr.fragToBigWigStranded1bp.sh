#!/usr/bin/env bash

# command-line application template

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [bed]
Description:
	Make a pair of stranded bigWig file from a fragment BED file using 1bp of each end
	- left-end => plus
	- right-end => minus
	in RPM scale (default) or manually scaled
Options:
	-o <outPrefix>: output prefix. required
		<outPrefix>.plus.bw
		<outPrefix>.minus.bw
	-g <chromSize>: chromosome size file, required
	-c <chrRegex>: Regular expression for chromosome selection, default=^chr[0-9XY]+$
		For multiple patterns use regular expression, such as \"^chr[0-9XY]+$|^chrM$\"
		NULL if not applicable or no filtering
	-m <memory>: memory size for sorting bedGraph file, default=5G
	-s <scale factor>: Manual scaling factor. This value is multiplied to \"raw read count\" primarily for spike-in based scaling.
			If 1, no scaling, i.e., raw read count.
			If 0, RPM normalized. default=0
	-n : Non-negative tracks
		In default, minus bigwig file contains negative values.
		But if this is set, it is generated in positive values like plus bigwig file.">&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
outPrefix=NULL
genome=NULL
chrRegex="^chr[0-9XY]+$"
memory=5G
nonnegative=FALSE
scaleFactor=0
while getopts ":o:g:c:m:s:n" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		c)
			chrRegex=$OPTARG
			;;
		m)
			memory=$OPTARG
			;;
		s)
			scaleFactor=$OPTARG
			;;
		n)
			nonnegative=TRUE
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

if [ "$genome" = "NULL" ];then
	echo -e "Error: genome (-g) must be specified" >&2
	exit 1
fi
assertFileExist $genome

if [ "$outPrefix" == "NULL" ];then
	echo -e "Error: Destination file must be specified" >&2
	exit 1
fi
desPlus=${outPrefix}.plus.bw
desMinus=${outPrefix}.minus.bw
if [ -f $desPlus ];then
	echo -e "Error: $desPlus already exists"
	exit 1
fi
if [ -f $desMinus ];then
	echo -e "Error: $desMinus already exists"
	exit 1
fi

src=$1
assertFileExist $src



###################################
## main code
printBed(){
	local src=$1
	assertFileExist $src

	local ext=${src##*.}
	
	if [ "$chrRegex" == "NULL" ];then
		if [ "$ext" == "gz" ];then
			zcat $src
		else
			cat $src
		fi
	else
		if [ "$ext" == "gz" ];then
			zcat $src | gawk '$1 ~ /'$chrRegex'/'
		else
			cat $src | gawk '$1 ~ /'$chrRegex'/'
		fi
	fi
}


desDir=`dirname $outPrefix`
mkdir -p $desDir

tmpBG=${TMPDIR}/__temp__.$$.bedGraph
tmpBwPlus=${TMPDIR}/__temp__.$$.plus.bw
tmpBwMinus=${TMPDIR}/__temp__.$$.minus.bw

echo -e "=============================================" >&2
echo -e "Creating BigWig file from a fragment bed file" >&2
echo -e "  - src = $src" >&2
echo -e "  - outPrefix = $outPrefix" >&2
echo -e "  - chromSize = $genome" >&2
echo -e "  - chrRegex = $chrRegex" >&2
echo -e "  - scaleFactor = $scaleFactor" >&2
echo -e "  - nonnegative = $nonnegative" >&2
echo -e "  - memory = $memory" >&2

if [ $scaleFactor == "0" ];then
	echo -e "  1) Calculating scale factor for RPM normalization" >&2
	ttc=`printBed $src | wc -l`
	scaleFactor=`echo $ttc | gawk '{ printf "%f", 1000000/$1}'`
	echo -e "\tTTC = $ttc (scaleFactor $scaleFactor)" >&2
else
	echo -e "  1) Scale factor was manually assigned" >&2
	echo -e "\tscaleFactor: ${scaleFactor}" >&2
fi

echo -e "  2) Making bedGraph file" >&2

echo -e "    - Plus strand" >&2
printBed $src \
	| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t+\n", $1,$2,$2+1,$4,$5 }' \
	| sort -S $memory -k1,1 -k2,2n -k3,3n \
	| genomeCoverageBed -bg -scale $scaleFactor -g $genome -i stdin \
	| gawk '{ printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4 }' \
	> $tmpBG
bedGraphToBigWig ${tmpBG} $genome ${tmpBwPlus}
rm $tmpBG

echo -e "    - Minus strand" >&2
if [ "$nonnegative" == "TRUE" ];then
	printBed $src \
		| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t+\n", $1,$3-1,$3,$4,$5 }' \
		| sort -S $memory -k1,1 -k2,2n -k3,3n \
		| genomeCoverageBed -bg -scale $scaleFactor -g $genome -i stdin \
		| gawk '{ printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4 }' \
		> $tmpBG
else
	printBed $src \
		| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t+\n", $1,$3-1,$3,$4,$5 }' \
		| sort -S $memory -k1,1 -k2,2n -k3,3n \
		| genomeCoverageBed -bg -scale $scaleFactor -g $genome -i stdin \
		| gawk '{ printf "%s\t%s\t%s\t-%.5f\n", $1,$2,$3,$4 }' \
		> $tmpBG
fi
bedGraphToBigWig ${tmpBG} $genome ${tmpBwMinus}
rm $tmpBG

mv ${tmpBwPlus} $desPlus
mv ${tmpBwMinus} $desMinus

