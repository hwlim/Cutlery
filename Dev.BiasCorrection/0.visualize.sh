#!/usr/bin/env bash

peak=/Volumes/limlab/Park/CnR.HNF4a/Default/3.Homer/Peak/JP1010_CnR_Hnf4a.homer.200bp.exBL.1rpm.bed

intersectBed -a ~/Research.office/Common_Data/MotifGenomeWide/Homer/mm10/HNF4a.bed -b $peak -u \
	| gawk '{ printf "%s\t%s\t%s\tHNF4a_%d\t%s\t%s\n", $1,$2,$3,NR,$5,$6 }' \
	> HNF4a_motif.in_peak.bed

marL=( 300 50 30 )
for mar in ${marL[@]}
do
	idom.visualizeExoBed.r -g homer_mm10 -f -t HNF4a_CnR -m 20,${mar} -s \
		-y 0.018,0.1 \
		-o HNF4a_motif.in_peak.${mar}bp \
		HNF4a_motif.in_peak.bed \
		test.corrected
done
#idom.visualizeExoBed.r -g mm10 -f -t HNF4a_CnR -m 20,300 -s \
#	-o HNF4a_motif.in_peak.300bp \
#	HNF4a_motif.in_peak.bed \
#	/Volumes/LimHeeWoong/Park/2.BigWig.stranded/JP1010_CnR_Hnf4a.filtered.dedup.nfr.sep
#
#idom.visualizeExoBed.r -g mm10 -f -t HNF4a_CnR -m 20,50 -s \
#	-o HNF4a_motif.in_peak.50bp \
#	HNF4a_motif.in_peak.bed \
#	/Volumes/LimHeeWoong/Park/2.BigWig.stranded/JP1010_CnR_Hnf4a.filtered.dedup.nfr.sep
#
#idom.visualizeExoBed.r -g mm10 -f -t HNF4a_CnR -m 20,30 -s \
#	-o HNF4a_motif.in_peak.30bp \
#	HNF4a_motif.in_peak.bed \
#	/Volumes/LimHeeWoong/Park/2.BigWig.stranded/JP1010_CnR_Hnf4a.filtered.dedup.nfr.sep
