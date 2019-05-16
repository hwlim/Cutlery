
########################
## Sample Information
import pandas as pd
samples = pd.read_csv("sample.tsv", sep="\t", comment="#")
#samples_indexById = samples.set_index("Id")
#samples_indexByName = samples.set_index("Name")

########################
genome="mm10"
chrom_size="chrom.size"
adapter="AGATCGGAAGAGC"
#trim_maxLen=100
trim_minLen=25
trim_minQual=20

star_index="index dir"


## Directories
fastqDir = "0.Fastq"
trimDir = fastqDir + "/Trim"
alignDir = "1.1.Align"
filteredDir = "1.2.Align.filtered"
dedupDir = "1.3.Align.dedup"
splitDir = "1.4.Align.split"
baseFreqDir = filteredDir + "/BaseFreq"
bigWigDir = "2.BigWig"

#def getfq(wildcards):
#	return "0.Fastq/" + samples_indexById.loc[wildcards.sampleId, ["Fq1","Fq2"]]
#
#def getfq_trim(wildcards):
#	sampleId = samples_indexByName.loc[wildcards.sampleName, ["Id"]]
#	return [ "0.Fastq/Trim/" + sampleId + "_1.trim.fq.gz", "0.Fastq/Trim/" + sampleId + "_2.trim.fq.gz" ]

	
rule all:
	input:
		expand(filteredDir + "/BaseFreq/{sampleName}.R{read}.freq.line.png",
			sampleName=samples.Name.tolist(), read=[1,2]),
		expand(bigWigDir + "/{sampleName}.filtered.dedup.{group}.ctr.bw",
			sampleName=samples.Name.tolist(), group=["nfr","nuc"])
#		expand("1.4.Align.split/{sampleId}.filtered.dedup.{group}.{proctype}.bed", 
#			sampleId=samples["Id"].tolist(),
#			group=["nfr","nuc"],
#			proctype=["con","ctr","sep"])
#		expand("1.3.Align.dedup/{sampleId}.filtered.dedup.bam", sampleId=samples["Id"].tolist()),
#		expand("1.2.Align.filtered/{sampleId}.filtered.bam", sampleId=samples.index.values.tolist()),

rule clean:
	shell:
		"rm -rf " + " ".join([ trimDir, alignDir, filteredDir, dedupDir, splitDir, baseFreqDir, bigWigDir ])

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
		minQual = trim_minQual,
	log:
		trimDir + "/{sampleId}.trim.log"
	shell:
		# Needs to be implemented as a quantum transaction
		"cutadapt -a {params.adapter} -A {params.adapter} --minimum-length {params.minLen} -q {params.minQual}"
		" -o __temp__.$$.1.fq.gz -p __temp__.$$.2.fq.gz {input.fq1} {input.fq2} 2> {log};"
		"mv __temp__.$$.1.fq.gz {output.fq1};" 
		"mv __temp__.$$.2.fq.gz {output.fq2}" 

rule align_pe:
	input:
		fq1 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_1.trim.fq.gz",
		fq2 = lambda wildcards: trimDir + "/" + samples.Id[samples.Name == wildcards.sampleName] + "_2.trim.fq.gz"
	output:
		alignDir+"/{sampleName}.bam"
	message:
		"Aligning... [{wildcards.sampleName}]"
	params:
		index=star_index,
		otherOpt=""
	log:
		alignDir + "/{sampleName}.star.log"
	threads: 4
	shell:
		#"--limitBAMsortRAM 10000000000
		#"--genomeLoad LoadAndKeep
		"STAR --runMode alignReads --genomeDir {params.index} "
		"--genomeLoad NoSharedMemory "
		"--readFilesIn <( zcat {input.fq1} ) <( zcat {input.fq2} )"
		"--runThreadN {threads} "
		"--outSAMtype BAM Unsorted "
		"--outFilterMultimapNmax 1 "
		"--alignMatesGapMax 2000 "
		"--alignIntronMax 1 "
		"--outReadsUnmapped None "
		"--outFileNamePrefix __temp__.$$ 2> {log}; "
		"mv __temp__.$$Aligned.out.bam {output}"

rule filter_align:
	input:
		alignDir+"/{sampleName}.bam"
	output:
		filteredDir + "/{sampleName}.filtered.bam"
	message:
		"Filtering... [{wildcards.sampleName}]"
	shell:
		"touch {output}"

rule dedup_align:
	input:
		filteredDir + "/{sampleName}.filtered.bam"
	output:
		dedupDir + "/{sampleName}.filtered.dedup.bam"
	message:
		"Deduplicating... [{wildcards.sampleName}]"
	shell:
		"touch {output}"

rule check_baseFreq:
	input:
		filteredDir + "/{sampleName}.filtered.bam"
	output:
		expand(baseFreqDir + "/{{sampleName}}.{read}.freq.line.png",
			read=["R1","R2"])
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"touch {output}"

rule split_bam:
	input:
		dedupDir + "/{sampleName}.filtered.dedup.bam"
	output:
		expand(splitDir + "/{{sampleName}}.filtered.dedup.{group}.{proctype}.bed",
			group=["nfr","nuc"], proctype=["con","ctr","sep"])
	message:
		"Splitting BAM file by fragment size... [{wildcards.sampleName}]"
	shell:
		"touch {output}"

rule make_bigwig:
	input:
		expand(splitDir + "/{{sampleName}}.filtered.dedup.{group}.ctr.bed",
			group=["nfr","nuc"])
	output:
		expand(bigWigDir + "/{{sampleName}}.filtered.dedup.{group}.ctr.bw",
			group=["nfr","nuc"])
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"touch {output}"

#rule peak_call:
#	input:
#	output:
#	message:
#	shell:
#
