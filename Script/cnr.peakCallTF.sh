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
	- <outDir>/peak.txt              Homer peak calling result
	- <outDir>/peak.bed              Homer peak in bed format
	- <outDir>/peak.exBL.bed         After blacklist filtering
	- <outDir>/peak.exBL.1rpm.bed    > 1rpm after filtering
Options:
	-o <outDir>: Destination tag directory, required
	-i <ctrl>: (optional) ctrl homer tag directory, default=NULL
	-m <mask>: mask bed file for filtering such as ENCODE blacklist, default=NULL
	-b <bigwig>: bigwig file of NFR fragments, i.e. *.nfr.con.bw
		If specified, peak centering is performed using given bigwig file
		default=NULL
	-s <optStr>: additional option for 'findPeaks' of Homer
		Internal pre-set option: \"-style factor -tbp 0 -inputtbp 0 -norm 1000000 -strand both -center -size 200 -C 0\"
		\"-fragLength 100 -inputFragLength 100\" is recommended as an additiona option considering the input bed length" >&2
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
bw=NULL
optStr=""
#lengthParam=NULL,NULL
while getopts ":o:i:m:b:s:" opt; do
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
		b)
			bw=$OPTARG
			;;
		s)
			optStr=$OPTARG
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

if [ "$bw" != "NULL" ];then
	assertFileExist $bw
fi

###################################
## main code
log=${desDir}/peak.log
echo -e "Homer peak-calling" >&2
echo -e "- target = $target" >&2
echo -e "- ctrl = $ctrl" >&2
echo -e "- desDir = $desDir" >&2
echo -e "- optStr = $optStr" >&2
echo -e "" >&2

peak0=${desDir}/peak.txt
peakBed=${desDir}/peak.bed
peakMasked=${desDir}/peak.exBL.bed
peak1rpm=${desDir}/peak.exBL.1rpm.bed
peakStat=${desDir}/peak.exBL.1rpm.stat

## Homer Peak finding
mkdir -p $desDir
optStr="-o ${peak0} -style factor -tbp 0 -inputtbp 0 -norm 1000000 -strand both -center -size 200 -C 0 ${optStr}"
if [ "$ctrl" == "NULL" ];then
	## Without control sample
	findPeaks $target ${optStr} 2>&1 | tee ${log}
else
	## With control sample
	findPeaks $target -i ${ctrl} ${optStr} 2>&1 | tee ${log}
fi

tmpAll=${TMPDIR}/__temp__.$$.all.bed
tmpCtr=${TMPDIR}/__temp__.$$.ctr.bed
tmpStat=${TMPDIR}/__temp__.$$.stat

## pos -> bed file
## - Filtering by $2 >1 : homer peak calling gives start coordinate 1 at the chromosome starting boundary with different size, which should be removed
grep -v "^#" ${peak0} \
	| gawk '{ printf "%s\t%d\t%d\t%s\t%s\t+\n", $2, $3, $4, $1, $6 }' \
	| gawk '$2 > 1' \
	> ${tmpAll}

## Peak centering using given NFR bigwig file
if [ "$bw" != "NULL" ];then
	cnr.centerPeaks.r -n 20 -o ${tmpCtr} ${tmpAll} ${bw}
	cat ${tmpCtr} | gawk '$2 > 0' > ${peakBed}
	rm $tmpCtr
else
	cat ${tmpAll} > ${peakBed}
fi
rm $tmpAll


## Blacklist filtering if given
if [ "$mask" != "NULL" ];then
	intersectBed -a ${peakBed} -b $mask -v > ${peakMasked}
else
	cat ${peakBed} > ${peakMasked}
fi

## > 1 RPM filtering
gawk '$5 > 1' ${peakMasked} > ${peak1rpm}

## Final summary print
N0=`cat ${peakBed} | wc -l`
N1=`cat ${peak1rpm} | wc -l`
echo -e "Final peaks: $N1 (/$N0), > 1rpm(/all)" >&2

## peak statistics, % of fragments in peak > 1rpm
N_frag_all=`grep genome ${target}/tagInfo.txt | gawk '{ printf "%d", $3}'`
tagDir2bed.pl $target -separate \
	| intersectBed -a stdin -b ${peak1rpm} -u \
	| wc -l \
	| gawk '{ printf "FragmentInPeakFrac\t%.2f\n", $1 / '$N_frag_all' * 100 }' \
	> $tmpStat
mv $tmpStat $peakStat
