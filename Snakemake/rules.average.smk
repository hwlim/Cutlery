########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
########$$$$$$$###################################


#####################################################
## Rules for averaging by "Group" column
def get_bigwig_rep(groupName, fragment):
	repL = samples.Name[samples.Group == groupName].tolist()
	return map(lambda x: bigWigDir + "/" + x + "." + fragment ".ctr.bw", repL)


rule make_bigwig_avg:
	input:
		all = lambda wildcards: get_bigwig_rep(wildcards.groupName, "all"),
		nfr = lambda wildcards: get_bigwig_rep(wildcards.groupName, "nfr"),
		nuc = lambda wildcards: get_bigwig_rep(wildcards.groupName, "nuc")
#		nfr = get_bigwig_rep_nfr,
#		nuc = get_bigwig_rep_nuc
	output:
		all = bigWigDir_avg + "/{groupName}.all.ctr.bw",
		nfr = bigWigDir_avg + "/{groupName}.nfr.ctr.bw",
		nuc = bigWigDir_avg + "/{groupName}.nuc.ctr.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		makeBigWigAverage.sh -g {chrom_size} -m 5G -o {output.all} {input.all}
		makeBigWigAverage.sh -g {chrom_size} -m 5G -o {output.nfr} {input.nfr}
		makeBigWigAverage.sh -g {chrom_size} -m 5G -o {output.nuc} {input.nuc}
		"""

