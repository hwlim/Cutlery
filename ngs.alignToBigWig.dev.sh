#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

# tagDir:	
# Name:		Data name or annotation
# Genome:	Genome
# fragLen:	Estimated fragment length
configL=(
3.Homer/JP1009_CnR_Ctrl.filtered.dedup.nfr.sep/TSV	JP1009_CnR_Ctrl.filtered.dedup.nfr.sep	mm10	1
3.Homer/JP1010_CnR_Hnf4a.filtered.dedup.nfr.sep/TSV	JP1010_CnR_Hnf4a.filtered.dedup.nfr.sep	mm10	1
)

src=<src>
desPrefix=<desPrefix>
chrom=<chrom>
fragLen=<fragLen>
sortMem=<sortMem>


assertFileExist $chrom

desDir=`dirname $desPrefix`
mkdir -p $desDir



printAlign(){
	if [ -d "$src" ];then
		## Homer TagDir
		echo -e "$src is a directory; assume it as a homer tag directory" >&2
		tagDir2bed.pl -separate $src
	elif [ -f "$src" ];then
		## BED(.gz) or BAM file
		ext=${src##*.}
		if [ "$ext" == "gz" ];then
			zcat $src
		elif [ "$ext" == "bed" ];then
			cat $src
		elif [ "$ext" == "bam" ];then
			bamToBed -i $src
		else
			echo -e "Error: $src is unknown format" >&2
			exit 1
		fi
	else
		echo -e "Error: $src does not exist" >&2
		exit 1
	fi
}



bed=${TMPDIR}/$$.bed.gz
tmpBg_plus=${TMPDIR}/$$.plus.bedGraph
tmpBg_minus=${TMPDIR}/$$.minus.bedGraph
tmpBw_plus=${TMPDIR}/$$.plus.bw
tmpBw_minus=${TMPDIR}/$$.minus.bw
bw_plus=${name}.plus.bw
bw_minus=${name}.minus.bw


## Resized alignment bed file
echo -e "1) Making (resized) bed file: fragLen=${fragLen}" >&2
if [ ${fragLen} -gt 0 ];then
	printAlign ${src} \
		| gawk '{
				if( $6 == "+" ){
					s = $2
					e = s + '$fragLen'
				}else{
					e = $3
					s = e - '$fragLen'
				}
				printf "%s\t%d\t%d\t.\t0\t%s\n", $1, s, e, $6
			}' \
		| sort -S $sortMem -k1,1 -k2,2n -k3,3n \
		> $bed
else
	printAlign ${src} \
		| sort -S $sortMem -k1,1 -k2,2n -k3,3n \
		> $bed
fi


## Total tag count
ttc=`zcat $bed | wc -l`

## BedGraph file
echo -e "2) Making bedGraph files\n\t TotalTags: $ttc" >&2
echo -e "\tStep 1: Generating Plus: $tmpBg_plus" >&2
genomeCoverageBed -bg -i ${bed} -strand + -g $chrom \
	| gawk 'BEGIN{ttc='$ttc'}
		{
			printf "%s\t%s\t%s\t%.5f\n", $1,$2,$3,$4*1000000/ttc;
		}'\
	> ${tmpBg_plus}
echo -e "\tStep 2: Generating Minus: $tmpBg_minus" >&2
genomeCoverageBed -bg -i ${bed} -strand - -g $chrom \
	| gawk 'BEGIN{ttc='$ttc'}
		{
			printf "%s\t%s\t%s\t-%.5f\n", $1,$2,$3,$4*1000000/ttc;
		}'\
	> ${tmpBg_minus}


## BigWig file
isFileExist $tmpBw_plus $tmpBw_minus
echo -e "3) Making bigWig files" >&2
echo -e "\tStep 1: Generating Plus: $bw_plus" >&2
bedGraphToBigWig ${tmpBg_plus} ${chrom} ${tmpBw_plus}
echo -e "\tStep 2: Generating Minus: $bw_minus" >&2
bedGraphToBigWig ${tmpBg_minus} ${chrom} ${tmpBw_minus}

mv ${tmpBw_plus} ${bw_plus}
mv ${tmpBw_minus} ${bw_minus}
