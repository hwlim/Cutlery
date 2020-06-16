########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
##################################################
'''
Required Variables
- adapter
- trim_minLen / trim_minQual
- star_index / star_option
- fastqDir
- trimDir
- alignDir
- filteredDir
- dedupDir
- sample [ Id, Name ]
'''

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
		#fq1 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_1.trim.fq.gz",
		#fq2 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_2.trim.fq.gz"
	output:
		alignDir+"/{sampleName}/align.bam"
	message:
		"Aligning... [{wildcards.sampleName}]"
	params:
		index=star_index,
		option=star_option
	log:
		alignDir + "/{sampleName}/star.log"
	threads:
		cluster["align_pe"]["cpu"]
	shell:
		"""
		module load CnR/1.0
		star.align.sh -g {params.index} \
			-o {alignDir}/{wildcards.sampleName}/align \
			-t {threads} \
			-p '{params.option}' \
			{input}
		"""
#		STAR --runMode alignReads --genomeDir {params.index} \
#			--genomeLoad NoSharedMemory \
#			--readFilesIn <( zcat {input.fq1} ) <( zcat {input.fq2} ) \
#			--runThreadN {threads} \
#			{params.option} \
#			--outFileNamePrefix __temp__.$$ 2> {log}
#		mv __temp__.$$Aligned.out.bam {output}"

rule filter_align:
	input:
		alignDir+"/{sampleName}/align.bam"
	output:
		filteredDir + "/{sampleName}.bam"
	message:
		"Filtering... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.filterBam.sh  -o {output} -c "{chrRegexAll}" {input}
		"""

rule dedup_align:
	input:
		filteredDir + "/{sampleName}.bam"
	output:
		dedupDir + "/{sampleName}.bam"
	message:
		"Deduplicating... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["dedup_align"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		cnr.dedupBam.sh -m {params.memory} -o {output} -r {input}
		"""
