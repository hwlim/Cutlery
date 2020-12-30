########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### for Processing pooled replicates
###
### - except for "rule all"
### 
###		Written by Hee Woong Lim
########$$$$$$$###################################

## Replicate-pooling -> bam file
def get_bam_replicate(wildcards):
	repL = sampleAll.Name[sampleAll.Group == wildcards.groupName].tolist()
	return map(lambda x: bamDir_rep + "/" + x + ".bam", repL)

rule pool_replicate_bam:
	input:
		get_bam_replicate
	output:
		bamDir_pool + "/{groupName}.bam"
	message:
		"Pooling replicates... [{wildcards.groupName}]"
	shell:
		"""
		module load Cutlery/1.0
		ngs.concateBamFiles.sh -o {output} {input}
		"""
