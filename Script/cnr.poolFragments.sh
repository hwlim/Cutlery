#!/usr/bin/env bash

###########################################3
# Cut&Run tools
# Written by Hee-Wooong Lim
# 
# Pool fragments files of multiple replicates by group

source $COMMON_LIB_BASE/commonBash.sh

function printUsage {
	echo -e "Usage: `basename $0` (options) <sample.tsv> <src bam directory> <des bam directory>
Description:
	Merge multiple fragment bed files of a group according to sample/group information within a given sample.tsv file
Input:
	- sample.tsv file: containing columns 'Name' and 'Group', from original/replicate Cutlery run
	- src replicate sample directory
	  e.g. <sample dir>
	   <sample dir>/sample1
	  ├── sample1
	  │    └── fragment.bed.gz
	  └── sample12
	       └── fragment.bed.gz
	- des directory pooling replicate samples per group
	  e.g. <des dir>
	  ├── group1
	  │    └── fragment.bed.gz
	  └── group2
	       └── fragment.bed.gz
Options:
	-n: number of cpu to use, default=1
	-m: memory to use for bsub, default=10000
	-t: if set, dry run simply displaying pooling message, default=off
	-b: if set, bsub are submitted for merging bam files, default=off
	-f: if set, force overwrite existing destination bam files, default=off" >&2
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

tmpPrefix=${TMPDIR}/__temp__.pool_frag.$$.$RANDOM
trap 'if [ `ls -1 ${tmpPrefix}* 2>/dev/null | wc -l` -gt 0 ];then rm ${tmpPrefix}*; fi' EXIT


###################################
## option and input file handling
memory=10000
cpu=1
testOnly=FALSE
bsub=FALSE
unsorted=FALSE
overwrite=FALSE
while getopts ":m:n:bfut" opt; do
	case $opt in
		m)
			memory=$OPTARG
			;;
		n)
			cpu=$OPTARG
			;;
		t)
			testOnly=TRUE
			;;
		b)
			bsub=TRUE
			;;
		u)
			unsorted=TRUE
			;;
		f)
			overwrite=TRUE
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
if [ $# -lt 3 ];then
	printUsage
	exit 1
fi
sampleInfo=$1
srcDir=$2
desDir=$3

assertFileExist $sampleInfo
assertDirExist $srcDir

###################################
## main code

groupL=`tail -n +2 $sampleInfo | grep -v -e ^$ -e "^#" -e ^Id | cut -f 3 | sort | uniq`

echo -e "Pooling replicate fragment files" >&2
echo -e "  - sampleInfo: $sampleInfo" >&2
echo -e "  - srcDir:   $srcDir" >&2
echo -e "  - desDir:   $desDir" >&2
echo -e "  - bsub:     $bsub" >&2
echo -e "  - unsorted: $unsorted" >&2
echo -e "  - testOnly: $testOnly" >&2

mkdir -p $desDir
for group in ${groupL[@]}
do
	des=${desDir}/${group}/fragment.bed.gz
	log=${desDir}/${group}/fragment.log

	## Checking existing destination file
	if [ -f $des ] && [ "$overwrite" != "TRUE" ];then
		echo -e "Warning: $des already exists, pass" >&2
		continue
	fi

	## List of replicate bam files
	srcL=( `cat $sampleInfo | grep -v -e ^$ -e "^#" -e ^Id | gawk '{ if($3 == "'$group'") printf "'$srcDir'/%s/fragment.bed.gz\n", $2 }'` )
	assertFileExist ${srcL[@]}

	## Merging to destination
	echo -e "Creating $des" >&2
	for src in ${srcL[@]}
	do
		echo -e "  - $src" >&2
	done

	## test only
	[ "$testOnly" == "TRUE" ] && continue


	## actual pooling
	mkdir -p ${desDir}/${group}
	echo -ne "" > $log
	for src in ${srcL[@]}
	do
		echo -e "- $src" >> $log
	done


	tmp=${tmpPrefix}_${group}.bed.gz
	if [ "$unsorted" == "TRUE" ];then
		if [ "$bsub" == "TRUE" ];then
			## Parallel processing using HPC:lsf
			bsub -W 24:00 -n 1 "cat ${srcL[@]} > $tmp; mv $tmp $des"
		else
			## Sequential processing
			cat ${srcL[@]} > $tmp
			mv $tmp $des
		fi
	else
		inputStr=""
		for src in ${srcL[@]}
		do
			inputStr="${inputStr} <( zcat $src )"
		done

		if [ "$bsub" == "TRUE" ];then
			## Parallel processing using HPC:lsf
			bsub -W 24:00 -n $cpu -M $memory -q rhel9 -R "span[hosts=1]" <<- EOF
#!/usr/bin/env bash
sort -m -k1,1 -k2,2n -k3,3n ${inputStr} | gzip > $tmp
mv $tmp $des
EOF
		else
			## Sequential processing
			#echo -e "sort -m -k1,1 -k2,2n -k3,3n ${inputStr}"
			eval "sort -m -k1,1 -k2,2n -k3,3n ${inputStr} | gzip > $tmp"
			mv $tmp $des
		fi
	fi
done


