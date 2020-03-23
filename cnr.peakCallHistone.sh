#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 
# Wrapper script for peak calling
#
# To do & consider:
# - default option & additional option handling

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) [taget tagDir]
Description: Make Homer data directory from BED file
Output:
	- <outDir>/peak.homer.txt              Homer peak calling result
	- <outDir>/peak.homer.bed              Homer peak in bed format
	- <outDir>/peak.homer.exBL.bed         After blacklist filtering
Options:
	-o <outDir>: Destination tag directory, required
	-i <ctrl>: (optional) ctrl homer tag directory, default=NULL
	-m <mask>: mask bed file for filtering such as ENCODE blacklist
	-s <optStr>: additional option for 'findPeaks' of Homer
		to internally pre-set option: \"-style histone -tbp 0 -norm 1000000 -strand both\"
		such as -size or -minDist
		\"-fragLength 100\" is recommended as an additiona option considering the input bed length" >&2
#	- <outDir>/peak.homer.exBL.1rpm.bed    > 1rpm after filtering
#	echo -e "\t-l <fragLen,peakWidth>: Comma-separated fragment length and peak width, default=NULL,NULL" >&2
#        echo -e "\t-g <genome>: genome, default=NULL" >&2
#        echo -e "\t-h: Print help" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
desDir=NULL
ctrl=NULL
mask=NULL
optStr="-style histone -tbp 0 -norm 1000000 -strand both"
#lengthParam=NULL,NULL
while getopts ":o:i:m:s:" opt; do
	case $opt in
		o)
			desDir=$OPTARG
			;;
		i)
			ctrl=$OPTARG
			;;
		m)
			mask=$OPTARG
			;;
		s)
			optStr="$optStr $OPTARG"
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

target=$1
assertDirExist $target

if [ "$desDir" == "NULL" ];then
	echo -e "Error: Destination directory (-o) must be specified" >&2
	exit 1
fi

if [ "$ctrl" != "NULL" ];then
	assertDirExist $ctrl
fi

if [ "$mask" != "NULL" ];then
	assertFileExist $mask
fi

assertFileExist ${target}/tagInfo.txt
ttc=`( grep genome ${target}/tagInfo.txt || true ) | cut -f 3`
if [ "$ttc" == "" ];then
	echo -e "Error: no tag count information found in [${target}/tagInfo.txt]" >&2
	exit 1
fi

###################################
## main code
log=${desDir}/peak.homer.log
echo -e "Homer peak-calling" >&2
echo -e "  - target = $target" >&2
echo -e "  - ctrl = $ctrl" >&2
echo -e "  - desDir = $desDir" >&2
echo -e "  - TTC = $ttc" >&2
echo -e "  - optStr = $optStr" >&2
echo -e "" >&2

peak0=${desDir}/peak.homer.txt
peakBed=${desDir}/peak.homer.bed
peakMasked=${desDir}/peak.homer.exBL.bed
#peak1rpm=${desDir}/peak.homer.exBL.1rpm.bed
tmpPeakMasked=${TMPDIR}/__temp__.$$.bed
tmpTagCount=${TMPDIR}/__temp__.$$.target

mkdir -p $desDir
if [ "$ctrl" == "NULL" ];then
	echo -e "findPeaks $target -o ${peak0} ${optStr}" >&2
	findPeaks $target -o ${peak0} ${optStr} 2>&1 | tee ${log}
else
	echo -e "findPeaks $target -i ${ctrl} -o ${peak0} ${optStr}" >&2
	findPeaks $target -i ${ctrl} -o ${peak0} ${optStr} 2>&1 | tee ${log}
fi


## To implement, how to calculate RPKM value efficiently for 5th column


echo -e "Convergting to bed" >&2
grep -v "^#" ${peak0} \
	| gawk '{ printf "%s\t%d\t%d\t%s\t%.3f\t+\n", $2, $3, $4, $1, $6*1000/($4-$3) }' \
	> ${peakBed}

echo -e "Blacking masking & merging" >&2
if [ "$mask" != "NULL" ];then
	subtractBed -a ${peakBed} -b $mask \
		| sortBed \
		| mergeBed \
		| gawk '{ printf "%s\t%d\t%d\tpeak.%d\t0\t+\n", $1,$2,$3,NR }' \
		> ${tmpPeakMasked}
		#| sort -k4,4 \
else
	cat ${peakBed} \
		| sortBed \
		| mergeBed \
		| gawk '{ printf "%s\t%d\t%d\tpeak.%d\t0\t+\n", $1,$2,$3,NR }' \
		> ${tmpPeakMasked}
fi

if [ `cat $tmpPeakMasked | wc -l` -eq 0 ];then
	echo -e "No remaining peaks" >&2
	touch $peakMasked
else
	echo -e "Tag counts in RPKM" >&2

	getPeakTags $tmpPeakMasked $target -tagAdjust 0 -tbp 0 -fixed \
		| sort -k1,1 \
		> ${tmpTagCount}

	paste ${tmpPeakMasked} ${tmpTagCount} \
		| gawk '{ printf "%s\t%d\t%d\tpeak.%d\t%.5f\t%s\n", $1,$2,$3,NR,$8*1000000/'${ttc}'*1000/($3-$2),$6 }' \
		| sort -k5,5nr \
		> $peakMasked
fi
cp ${tmpPeakMasked} .
cp ${tmpTagCount} .
#if [ "$ctrl" != "NULL" ];then
#	getPeakTags $peakMasked $ctrl -tagAdjust 0 -tbp 0 \
#		| sort -k1,1 \
#		> __temp__.$$.ctrl
#fi

#gawk '$5 > 1' ${peakMasked} > ${peak1rpm}

N0=`cat ${peakBed} | wc -l`
N1=`cat ${peakMasked} | wc -l`
echo -e "Final peaks:" >&2
echo -e "  - Original = $N0" >&2
echo -e "  - After filtering = $N1" >&2
