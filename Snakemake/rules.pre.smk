########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
##################################################


## default STAR module
if 'star_module' not in locals():
	star_module = "STAR/2.7.4"

if 'opt_cutadapt' not in locals():
	opt_cutadapt = ""

if 'chrRegexAll' not in locals():
	chrRegexAll = chrRegexTarget


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
		module purge
		module load cutadapt/2.1.0
		cutadapt -a {params.adapter} -A {params.adapter} --minimum-length {params.minLen} -q {params.minQual} {opt_cutadapt} \
			-o __temp__.$$.1.fq.gz -p __temp__.$$.2.fq.gz {input.fq1} {input.fq2} 2>&1 | tee {log}
		mv __temp__.$$.1.fq.gz {output.fq1}
		mv __temp__.$$.2.fq.gz {output.fq2} 
		"""

def get_fastq(sampleName):
	if doTrim:
		return [trimDir + "/" + samples.Id[samples.Name == sampleName].tolist()[0] + "_1.trim.fq.gz",
			trimDir + "/" + samples.Id[samples.Name == sampleName].tolist()[0] + "_2.trim.fq.gz"]
	else:
		fq1=samples.Fq1[samples.Name == sampleName].tolist()[0]
		fq2=samples.Fq2[samples.Name == sampleName].tolist()[0]
		if(fq2=="NULL"):
			return fastqDir + "/" + fq1
		else:
			return [fastqDir + "/" + fq1, fastqDir + "/" + fq2]

# def get_fq1(sampleName):
# 	if doTrim:
# 		return trimDir + "/" + samples.Id[samples.Name == sampleName].tolist()[0] + "_1.trim.fq.gz"
# 	else:
# 		return fastqDir + "/" + samples.Fq1[samples.Name == sampleName].tolist()[0]

# def get_fq2(sampleName):
# 	if doTrim:
# 		return trimDir + "/" + samples.Id[samples.Name == sampleName].tolist()[0] + "_2.trim.fq.gz"
# 	else:
# 		return fastqDir + "/" + samples.Fq2[samples.Name == sampleName].tolist()[0]


rule align_pe:
	input:
		fq = lambda wildcards: get_fastq(wildcards.sampleName)
	output:
		bam = alignDir + "/{sampleName}/align.sortByName.bam"
		#bai = alignDir + "/{sampleName}/align.bam.bai"
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
		module purge
		module load Cutlery/1.0
		module load {params.star_module}

		STAR --runMode alignReads \
			--genomeDir {params.index} \
			--readFilesIn {input.fq} \
			--readFilesCommand zcat \
			--genomeLoad NoSharedMemory \
			--outSAMtype BAM Unsorted \
			--outFileNamePrefix {alignDir}/{wildcards.sampleName}/align. \
			--runThreadN {threads} \
			--outTmpDir ${{TMPDIR}}/STARtmp_$$_$RANDOM \
			{star_option}
		
		#	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 \
		#mv {alignDir}/{wildcards.sampleName}/align.Aligned.sortedByCoord.out.bam {alignDir}/{wildcards.sampleName}/align.bam
		#samtools index {alignDir}/{wildcards.sampleName}/align.bam

		mv {alignDir}/{wildcards.sampleName}/align.Aligned.out.bam {alignDir}/{wildcards.sampleName}/align.sortByName.bam

		if [ -f {alignDir}/{wildcards.sampleName}/align.Unmapped.out.mate1 ];then
			gzip {alignDir}/{wildcards.sampleName}/align.Unmapped.out.mate1
			gzip {alignDir}/{wildcards.sampleName}/align.Unmapped.out.mate2
		fi
		"""


rule csort_bam:
	input:
		bam=alignDir+"/{sampleName}/align.sortByName.bam"
	output:
		bam = alignDir+"/{sampleName}/align.bam",
		bai = alignDir+"/{sampleName}/align.bam.bai"
	message:
		"Sorting by coordinate... [{wildcards.sampleName}]"
	threads:
		cluster["csort_bam"]["cpu"]
	shell:
		"""
		module purge
		module load ChIPseq/1.0
		samtools sort -o {output.bam} -T ${{TMPDIR}}/csort_bam.{wildcards.sampleName}.${{RANDOM}} -@ {threads} -m 2G {input.bam}
		samtools index {output.bam}
		"""


def get_align_dir(bamList):
	import os.path
	return list(map(lambda x: os.path.dirname(x), bamList ))

rule make_align_stat_table:
	input:
		expand(alignDir+"/{sampleName}/align.sortByName.bam", sampleName=samples.Name.tolist())
	output:
		qcDir + "/alignStat.txt"
	params:
		inputDir = get_align_dir(expand(alignDir+"/{sampleName}/align.sortByName.bam", sampleName=samples.Name.tolist())),
		outPrefix = lambda wildcards, output: __import__("re").sub(".txt$","", output[0])
	message:
		"Creating alignment stat file"
	shell:
		"""
		module purge
		module load R/4.4.0
		star.getAlignStats.r -o {params.outPrefix} {params.inputDir}
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


##########################################
## Rule to handle multimapper using CSEM

## Run CSEM for postprocessing of multimapper
# 20240710: better to make output temporary by Snakemake
rule run_csem:
	input:
		bam = alignDir+"/{sampleName}/align.sortByName.bam"
	output:
		bam = alignDir+"/{sampleName}/CSEM/align.bam"
	message:
		"Running CSEM... [{wildcards.sampleName}]"
	params:
		#outDir = lambda wildcards, output: __import__("os").path.dirname(output[0]),
		prefix = lambda wildcards, output: __import__("re").sub(".bam$", "", output[0])
	threads:
		cluster["run_csem"]["cpu"]
	shell:
		"""
		#module load ChIPseq/1.0
		module purge
		module load R/4.4.0
		module load bedtools/2.30.0
		module load csem_limlab/06272024
		ngs.run_csem.sh -o {params.prefix} -t {threads} -n {wildcards.sampleName} {input.bam}
		"""

## Unify CSEM bam file
rule unify_csem:
	input:
		alignDir+"/{sampleName}/CSEM/align.bam"
	output:
		bam = alignDir+"/{sampleName}/CSEM/align.uniq.bam",
		bai = alignDir+"/{sampleName}/CSEM/align.uniq.bam.bai"
	message:
		"Unifying CSEM results... [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0]),
		prefix = lambda wildcards, output: __import__("re").sub(".bam$", "", output[0])
	shell:
		"""
		module purge
		module load ChIPseq/1.0
		tmpBam=${{TMPDIR}}/unify_csem.{wildcards.sampleName}.${{RANDOM}}.bam

		ngs.unifyCSEM.py -o ${{tmpBam}} {input}
		samtools sort -o {output.bam} -T ${{TMPDIR}}/unify_csem.sort.{wildcards.sampleName}.${{RANDOM}} -@ {threads} -m 2G ${{tmpBam}}
		rm ${{tmpBam}}
		#mv ${{tmpBam}} {output.bam}
		samtools index {output.bam}
		"""

