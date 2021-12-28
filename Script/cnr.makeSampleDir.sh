#!/usr/bin/env bash

## Convert existing Cutlery processing results into the current version of sample folders
## Organizing all the output files in sample by sample manner not by output type 
## To save storage, existing out files/folders do not change. Instead, symbolic links are generated

source $COMMON_LIB_BASE/commonBash.sh


getDir()
{
	# $1 : Snakefile
	# $2 : variable name to retrieve
	src=$1
	name=$2

	set +o pipefail
	result=`head -n 1000 $src | grep ${name} | head -n 1 | sed -e 's/[ \t]*//g' -e 's/=/\t/' -e 's/"//g' | cut -f 2`
	set -o pipefail

	if [ "$result" == "" ];then
		echo -e "Error: no variable named [$name] in $src" >&2
		exit 1
	fi

	echo $result
}

makeLink()
{
	src=$1
	des=$2

	if [ ! -f $src ] && [ ! -d $src ];then
		echo -e "  Warning: $src does not exists, skip" >&2
		return
	fi

	if [ -f $des ] || [ -d $des ];then
		echo -e "  Warning: $des already exists, skip" >&2
	else
		ln -v -s $src $des
		#echo -e "  Creating a link: $des --> $src" >&2
	fi
}

nameL=(
bigWigDir1bp
bigWigDir1bp_abs
bigWigDirAllFrag
bigWigDir
bamDir_pool
homerDir
)

if [ $# -lt 2 ];then
	echo -e "Usage: `basename $0` <sample.tsv> <Snakefile>" >&2
	exit 1
fi

src_sample=$1
src_smk=$2
assertFileExist $src_sample
assertFileExist $src_smk

if [ ${src_sample} == "sample.tsv" ];then
	sampleL=( `tail -n +2 $src_sample | grep -v -e "^#" -e "^$" | cut -f 2` )
else
	sampleL=( `tail -n +2 $src_sample | grep -v -e "^#" -e "^$" | cut -f 1` )
fi

sampleDir=5.Samples
mkdir -p $sampleDir 

for sample in ${sampleL[@]}
do
	echo "Processing $sample" >&2
	desDir=${sampleDir}/${sample}
	mkdir -v -p $desDir

	echo -e "  Entering $desDir" >&2
	cd $desDir

	#for name in ${nameL[@]}
	#do
	#	dirName=`getDir ../../${src_smk} $name`
	#	echo $dirName >&2

	#bamDir_pool
	#splitDir
	#bigWigDir
	#bigWigDirAllFrag
	#bigWigDir1bp
	#bigWigDir1bp_abs
	#homerDir

	## bam
	dirName=`getDir ../../${src_smk} bamDir_pool`
	makeLink ../../${dirName}/${sample}.bam align.bam

	## fragments
	#mkdir -p Fragments
	#cd Fragments
	dirName=`getDir ../../${src_smk} splitDir`
	makeLink ../../${dirName}/${sample}.all.con.bed.gz fragment.bed.gz
	#makeLink ../../../${dirName}/${sample}.all.con.bed.gz frag.all.con.bed.gz
	#makeLink ../../../${dirName}/${sample}.nfr.con.bed.gz frag.nfr.con.bed.gz
	#makeLink ../../../${dirName}/${sample}.nuc.con.bed.gz frag.nuc.con.bed.gz
	#makeLink ../../../${dirName}/${sample}.all.ctr.bed.gz frag.all.ctr.bed.gz
	#makeLink ../../../${dirName}/${sample}.nfr.ctr.bed.gz frag.nfr.ctr.bed.gz
	#makeLink ../../../${dirName}/${sample}.nuc.ctr.bed.gz frag.nuc.ctr.bed.gz
	#makeLink ../../../${dirName}/${sample}.all.sep.bed.gz frag.all.sep.bed.gz
	#makeLink ../../../${dirName}/${sample}.nfr.sep.bed.gz frag.nfr.sep.bed.gz
	#makeLink ../../../${dirName}/${sample}.nuc.sep.bed.gz frag.nuc.sep.bed.gz
	#cd ..

	## Homer Tag Dir
	dirName=`getDir ../../${src_smk} homerDir`
	makeLink ../../${dirName}/${sample}/TSV.all TSV.all
	makeLink ../../${dirName}/${sample}/TSV.nuc TSV.nuc
	makeLink ../../${dirName}/${sample}/TSV.nfr TSV.nfr

	## bigwig ctr rpm
	dirName=`getDir ../../${src_smk} bigWigDir`
	makeLink ../../${dirName}/${sample}.all.ctr.bw igv.all.ctr.bw
	makeLink ../../${dirName}/${sample}.nfr.ctr.bw igv.nfr.ctr.bw
	makeLink ../../${dirName}/${sample}.nuc.ctr.bw igv.nuc.ctr.bw

	## bigwig allfrag rpm
	dirName=`getDir ../../${src_smk} bigWigDirAllFrag`
	makeLink ../../${dirName}/${sample}.allFrag.bw igv.all.con.bw
	makeLink ../../${dirName}/${sample}.nfr.con.bw igv.nfr.con.bw
	makeLink ../../${dirName}/${sample}.nuc.con.bw igv.nuc.con.bw

	## bigwig 1bp rpm
	dirName=`getDir ../../${src_smk} bigWigDir1bp`
	makeLink ../../${dirName}/${sample}.plus.bw igv.1bp.plus.bw
	makeLink ../../${dirName}/${sample}.minus.bw igv.1bp.minus.bw

	## bigwig 1bp raw/abs
	dirName=`getDir ../../${src_smk} bigWigDir1bp_abs`
	makeLink ../../${dirName}/${sample}.raw.abs.plus.bw igv.1bp.raw.abs.plus.bw
	makeLink ../../${dirName}/${sample}.raw.abs.minus.bw igv.1bp.raw.abs.minus.bw

	## peak
	dirName=`getDir ../../${src_smk} homerDir`
	makeLink ../../${dirName}/${sample}/HomerPeak.factor HomerPeak.factor
	makeLink ../../${dirName}/${sample}/HomerPeak.histone HomerPeak.histone

	echo -e "  Leaving $desDir" >&2
	cd ../../
done



