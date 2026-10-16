#!/usr/bin/env bash

###########################################
# Convert MACS2 --call-summits output (macs_peaks.xls) into summit-centered BED files
# with RPM-scaled summit pileup in column 5, mirroring Homer's peak.exBL.bed / peak.exBL.1rpm.bed

source $COMMON_LIB_BASE/commonBash.sh

function printUsage {
	echo -e "Usage: `basename $0` (options) [macs_peaks.xls]
Description: Make summit-centered peak BED files from MACS2 peak calling with --call-summits
Input:
	- macs_peaks.xls from 'macs2 callpeak --call-summits'
	  Requires the header line '# total fragments in treatment: <N>' (BAMPE mode)
Output:
	- <outPrefix>.exBL.bed         All summits after blacklist filtering; 5th column = summit pileup in RPM
	- <outPrefix>.exBL.1rpm.bed    > 1rpm after filtering
	Both are sorted by the 5th column in decreasing order
Options:
	-o <outPrefix>: Output prefix, required, e.g. MACS2.factor/macs_summits
	-m <mask>: mask bed file for filtering such as ENCODE blacklist, default=NULL
	-c <chrRegex>: regular expression of chromosomes to keep, default=\"^chr[0-9XY]+$\"
	-w <width>: peak width centered at the summit, default=200
		Summits closer than width/2 to the chromosome start are dropped
	-h: Print help" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
outPrefix=NULL
mask=NULL
chrRegex="^chr[0-9XY]+$"
width=200
while getopts ":o:m:c:w:h" opt; do
	case $opt in
		o)
			outPrefix=$OPTARG
			;;
		m)
			mask=$OPTARG
			;;
		c)
			chrRegex=$OPTARG
			;;
		w)
			width=$OPTARG
			;;
		h)
			printUsage
			exit 0
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
if [ $# -ne 1 ];then
	printUsage
	exit 1
fi
xls=$1

if [ "$outPrefix" == "NULL" ];then
	echo -e "Error: Output prefix (-o) must be specified" >&2
	exit 1
fi

assertFileExist $xls
if [ "$mask" != "NULL" ];then
	assertFileExist $mask
fi


###################################
## main code
tmpDir=`mktemp -d "${TMPDIR:-/tmp}/cnr.macsSummitToBed.XXXXXX"`
trap 'rm -rf "$tmpDir"' EXIT
set -eo pipefail

outAll=${outPrefix}.exBL.bed
out1rpm=${outPrefix}.exBL.1rpm.bed
tmpBed=${tmpDir}/summit.bed

echo -e "MACS summit to BED" >&2
echo -e "- xls = $xls" >&2
echo -e "- outPrefix = $outPrefix" >&2
echo -e "- mask = $mask" >&2
echo -e "- chrRegex = $chrRegex" >&2
echo -e "- width = $width" >&2
echo -e "" >&2

## Total fragments, e.g. "# total fragments in treatment: 6167439"
total=`gawk -F': ' '/^# total fragments in treatment:/ { print $2; exit }' $xls`
if [ -z "$total" ] || [ "$total" -eq 0 ];then
	echo -e "Error: '# total fragments in treatment' not found in $xls" >&2
	exit 1
fi

## xls columns: chr start end length abs_summit pileup -log10(p) fold -log10(q) name
## abs_summit is 1-based; 0-based summit = abs_summit - 1 (same position as macs_summits.bed)
grep -v '^#' $xls \
	| gawk -v total=$total -v width=$width -v re="$chrRegex" '
		BEGIN { hw = int(width / 2) }
		NF >= 10 && $1 != "chr" && $1 ~ re {
			s = $5 - 1
			if (s - hw < 0) next
			printf "%s\t%d\t%d\t%s\t%.4f\t+\n", $1, s - hw, s - hw + width, $10, $6 * 1000000 / total
		}' \
	> $tmpBed

## Blacklist filtering if given
if [ "$mask" != "NULL" ];then
	intersectBed -a $tmpBed -b $mask -v | sort -k5,5nr > $outAll
else
	sort -k5,5nr $tmpBed > $outAll
fi

## > 1 RPM filtering
gawk '$5 > 1' $outAll > $out1rpm

## Final summary print
N0=`cat $outAll | wc -l`
N1=`cat $out1rpm | wc -l`
echo -e "Final summits: $N1 (/$N0), > 1rpm(/all)" >&2
