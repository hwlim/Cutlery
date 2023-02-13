#!/usr/bin/env bash

###########################################3
# Written by Hee-Wooong Lim
# 

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo -e "Usage: `basename $0` (options) <src dir> <des dir>
Description:
	Create a sample folder containing a symbolic link to a fragment file
	To perform alternative Cutlery analysis from fragment file
Input:
	- src sample folder 
	- des sample folder
	
	A folder of sample folders that containing fragment files
	e.g.
	3.Sample/
	├── Ctrl
	├── dnRAR
	└── v2"
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi



###################################
## option and input file handling
#force=FALSE
#lengthParam=NULL,NULL
while getopts ":" opt; do
	case $opt in
#		f)
#			force=TRUE
#			;;
#		s)
#			optStr=$OPTARG
#			;;
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
if [ $# -lt 2 ];then
	printUsage
	exit 1
fi

srcDir="$1"
desDir="$2"

assertDirExist $srcDir
mkdir -p $desDir

for x in ${srcDir}/*/fragment.bed.gz
do
	name=`dirname $x`;
	name=`basename $name`;
	mkdir -p ${desDir}/${name};
	ln -srv ${x} ${desDir}/${name}/fragment.bed.gz;
done
