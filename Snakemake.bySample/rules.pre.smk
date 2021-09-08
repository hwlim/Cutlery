########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
##################################################


## default STAR module
if 'star_module' not in locals():
	star_module="STAR/2.5"


rule trim_pe:
	input:
		fq1 = lambda wildcards: fastqDir + "/" + samples.Fq1[samples.Id == wildcards.sampleId],
		fq2 = lambda wildcards: fastqDir + "/" + samples.Fq2[samples.Id == wildcards.sampleId]
	output:
		fq1 = trimDir + "/{sampleId}_1.trim.fq.gz",
		fq2 = trimDir + "/{sampleId}_2.trim.fq.gz"
	message:
		"Trimming... [{wildcards.sampleId}]"
	params:
		adapter = adapter,
		minLen = trim_minLen,
		minQual = trim_minQual
	log:
		trimDir + "/{sampleId}.trim.log"
	shell:
		"""
		# Note: Needs to be implemented as a quantum transaction
		cutadapt -a {params.adapter} -A {params.adapter} --minimum-length {params.minLen} -q {params.minQual} \
			-o __temp__.$$.1.fq.gz -p __temp__.$$.2.fq.gz {input.fq1} {input.fq2} 2>&1 | tee {log}
		mv __temp__.$$.1.fq.gz {output.fq1}
		mv __temp__.$$.2.fq.gz {output.fq2} 
		"""

def get_fastq(wildcards):
	if doTrim:
		return [trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName].tolist()[0] + "_1.trim.fq.gz",
			trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName].tolist()[0] + "_2.trim.fq.gz"]
	else:
		fq1=samples.Fq1[samples.Name == wildcards.sampleName].tolist()[0]
		fq2=samples.Fq2[samples.Name == wildcards.sampleName].tolist()[0]
		if(fq2=="NULL"):
			return fastqDir + "/" + fq1
		else:
			return [fastqDir + "/" + fq1, fastqDir + "/" + fq2]


rule align_pe:
	input:
		get_fastq
	output:
		bam = alignDir + "/{sampleName}/align.bam",
		bai = alignDir + "/{sampleName}/align.bam.bai"
	message:
		"Aligning... [{wildcards.sampleName}]"
	params:
		index = star_index,
		option = star_option,
		star_module = star_module
	log:
		alignDir + "/{sampleName}/star.log"
	threads:
		cluster["align_pe"]["cpu"]
	shell:
		"""
		module load Cutlery/1.0
		module load {params.star_module}

		star.align.sh -g {params.index} \
			-o {alignDir}/{wildcards.sampleName}/align \
			-t {threads} \
			-p '{params.option}' \
			-s \
			{input}
		"""

def get_align_dir(bamList):
	import os.path
	return list(map(lambda x: os.path.dirname(x), bamList ))

rule make_align_stat_table:
	input:
		expand(alignDir+"/{sampleName}/align.bam", sampleName=samples.Name.tolist())
	output:
		qcDir + "/alignStat.txt"
	params:
		inputDir = get_align_dir(expand(alignDir+"/{sampleName}/align.bam", sampleName=samples.Name.tolist()))
	message:
		"Creating alignment stat file"
	shell:
		"""
		module load ChIPseq/1.0
		star.getAlignStats.r {params.inputDir} > {output}
		"""

'''
## No more used since chromosome / flag filtering are now incorporated into bam to fragment conversion step
rule filter_align:
	input:
		alignDir + "/{sampleName}/align.bam"
	output:
		bam = filteredDir + "/{sampleName}.bam",
		bai = filteredDir + "/{sampleName}.bam.bai"
	message:
		"Filtering... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.filterBam.sh  -o {output.bam} -c "{chrRegexAll}" {input}
		samtools index {output.bam}
		"""
'''

rule dedup_align:
	input:
		alignDir + "/{sampleName}.bam"
	output:
		bam = dedupDir + "/{sampleName}.bam",
		bai = dedupDir + "/{sampleName}.bam.bai"
	message:
		"Deduplicating... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["dedup_align"]["memory"]/1000 - 1 )
	shell:
		"""
		module load Cutlery/1.0
		cnr.dedupBam.sh -m {params.memory} -o {output.bam} -r {input}
		samtools index {output.bam}
		"""
