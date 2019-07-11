#!/usr/bin/env bash

source $MYBASHLIB/commonBash.sh

srcFileL=(
SRR7965431_1.fastq.gz
SRR7965431_2.fastq.gz
SRR7965435_1.fastq.gz
SRR7965435_2.fastq.gz
)

n=4000

srcBase=/Volumes/limlab/PublicData/Collection.CnR_MNase/GSE111121_uliCnR/0.Fastq

set +o pipefail
for srcFile in ${srcFileL[@]}
do
	src=${srcBase}/${srcFile}
	des=${srcFile}
	if [ "$src" = "$des" ];then
		echo -e "Error: source and destination are the same" >&2
		exit 1
	fi

	echo -e "Creating $des" >&2
	zcat $src | head -n $n | gzip > $des
done
