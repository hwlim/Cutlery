#!/usr/bin/env bash
trap 'if [ `ls -1 __temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm __temp__.$$.*; fi' EXIT
source $MYBASHLIB/commonBash.sh


## k-mer configurations
getKmerEnrichment(){
	local src=$1
	local prefix=$2
	local kmer_up=$3
	local kmer_down=$4

	local kmer=`echo -e "${kmer_up}\t${kmer_down}" | gawk '{ print $1 + $2 }'`
	local kmerCntGenome
	if [ -d /home/heewlim/Research/Common_Data/${genome}/seqOutBias.autosome/ ];then
		kmerCntGenome=/home/heewlim/Research/Common_Data/${genome}/seqOutBias.autosome/kmer_count.read50.kmer${kmer}
	else
		kmerCntGenome=/home/heewlim/Research/Common_Data/${genome}/seqOutBias/kmer_count.read50.kmer${kmer}
	fi
	assertFileExist $kmerCntGenome


	local blist=/home/heewlim/Research/Common_Data/${genome}/ENCODE-blacklist.bed
	local desDir=.
	local kmerCnt=${desDir}/${prefix}.kmerCnt.txt
	local scaleFactor=${desDir}/${prefix}.scaleFactor.txt
	mkdir -p $desDir
	
	echo -e "Counting k-mers" >&2
	echo -e "  src : $src">&2
	echo -e "  prefix : $prefix">&2
	echo -e "  genome : $genome" >&2
	echo -e "  Kmer: ${kmer_up} / ${kmer_down}" >&2
	
	if [ -f ${kmerCnt} ];then
		echo -e "  $kmerCnt already exits, pass" >&2
	else
		echo -e "  Creating $kmerCnt" >&2
		zcat $src \
			| gawk '/^'$chrRx'/' \
			| intersectBed -a stdin -b $blist -v \
			| ~/Research/ExoTools/7.3.ExoKmer.selectChr/countKmers.sh -k ${kmer_up},${kmer_down} -g $genome stdin \
			> __temp__.$$.txt
		mv __temp__.$$.txt ${kmerCnt}
	fi
		
	if [ -f ${scaleFactor} ];then
		echo -e "  ${scaleFactor} already exists, pass" >&2
	else
		echo -e "  Creating $scaleFactor" >&2
		/home/heewlim/Research/ExoTools/7.3.ExoKmer.selectChr/calcKmerScaleFactor.sh ${kmerCntGenome} ${kmerCnt} > __temp__.$$.txt
		mv __temp__.$$.txt $scaleFactor
	fi
	echo -e "" >&2
}


startDir=`pwd`

for workDir in `ls -d */Up*/`
do
	echo -e "================================================" >&2
	echo -e "Entering $workDir" >&2
	cd $workDir
	source ../Background/config.sh
	source ./config.sh


	src=exo.fg.bed.gz

	if [ -f $exo ];then
		getKmerEnrichment $exo exo $kmer_up $kmer_down
	else
		echo -e "Warning $exo does not exists, pass" >&2
	fi

	if [ -f $chip ];then
		getKmerEnrichment $chip seq $kmer_up $kmer_down
	else
		echo -e "Warning $chip does not exists, pass" >&2
	fi

	cd $startDir
	echo -e "" >&2
done
	
