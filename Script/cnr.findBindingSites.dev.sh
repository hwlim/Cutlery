#!/usr/bin/env bash

srcPeak=Anchor/anchor.chr1.bed
srcFcl=hnf4a.cfl.bed.gz
window=400
desFclSelect=fclSelect.txt

extendBed.sh -s -w 300 $srcPeak \
	| gawk '{ printf "%s\t%d\t%d\tBED.%d\t%s\t%s\n", $1,$2,$3,NR,$5,$6 }' \
	| intersectBed -a stdin -b $srcFcl -wa -wb -sorted \
	| cut -f 1-4,8,11 \
	> $desFclSelect
