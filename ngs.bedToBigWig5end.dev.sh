#!/usr/bin/env bash

source $MYBASHLIB/commonBash.sh

# tagDir:	
# Name:		Data name or annotation
# Genome:	Genome
# fragLen:	Estimated fragment length
configL=(
3.Homer/JP1009_CnR_Ctrl.filtered.dedup.nfr.sep/TSV	JP1009_CnR_Ctrl.filtered.dedup.nfr.sep	mm10	1
3.Homer/JP1010_CnR_Hnf4a.filtered.dedup.nfr.sep/TSV	JP1010_CnR_Hnf4a.filtered.dedup.nfr.sep	mm10	1
)

srcBase=..
desDir=.
mkdir -p ${desDir}

for (( i=0;i<${#configL[@]};i=$i+4 ))
do
	tDir=${srcBase}/${configL[$i]}
	name=${configL[$i+1]}
	genome=${configL[$i+2]}
	fragLen=${configL[$i+3]}
	isDirExist $tDir

	bed=${name}.bed.gz
	bg_plus=${name}.plus.bedGraph
	bg_minus=${name}.minus.bedGraph
	bw_plus=${name}.plus.bw
	bw_minus=${name}.minus.bw


	echo -e "Generating bigWig files:" >&2
	echo -e "tDir=${tDir}\nData name=${name}\nGenome=${genome}\nFragLen=${fragLen}" >&2

	echo -e "1) Making extended bed file: fragLen=${fragLen}" >&2
	tagDir2bed.pl ${tDir} -len ${fragLen} -separate | gzip -c > $bed
	ttc=`zwcl.sh $bed | cut -f 1`

	isFileExist $bed
	echo -e "2) Making bedGraph files\n\t TotalTags: $ttc" >&2
	echo -e "\tStep 1: Generating Plus: $bg_plus" >&2
	genomeCoverageBed -bg -i ${bed} -strand + -g ~/Research/Common_Data/${genome}/chrom.sizes  \
		| gawk 'BEGIN{ttc='$ttc'}
			{
				printf "%s\t%s\t%s\t%.3f\n", $1,$2,$3,$4*1000000/ttc;
			}'\
		> ${bg_plus}
	echo -e "\tStep 2: Generating Minus: $bg_minus" >&2
	genomeCoverageBed -bg -i ${bed} -strand - -g ~/Research/Common_Data/${genome}/chrom.sizes  \
		| gawk 'BEGIN{ttc='$ttc'}
			{
				printf "%s\t%s\t%s\t-%.3f\n", $1,$2,$3,$4*1000000/ttc;
			}'\
		> ${bg_minus}


	isFileExist $bg_plus $bg_minus
	echo -e "3) Making bigWig files" >&2
	echo -e "\tStep 1: Generating Plus: $bw_plus" >&2
	bedGraphToBigWig ${bg_plus} ~/Research/Common_Data/${genome}/chrom.sizes ${bw_plus}
	echo -e "\tStep 2: Generating Minus: $bw_minus" >&2
	bedGraphToBigWig ${bg_minus} ~/Research/Common_Data/${genome}/chrom.sizes ${bw_minus}

#	echo -e "4) Compressing bedGraph files" >&2
#	gzip -v -c ${bg_plus} > ${bg_plus}.gz
#	gzip -v -c ${bg_minus} > ${bg_minus}.gz
	rm ${bg_plus} ${bg_minus}
done

