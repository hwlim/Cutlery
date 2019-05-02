import pandas as pd

samples = pd.read_csv("sample.tsv", sep="\t", comment="#").set_index("Id")

fastqDir="0.Fastq"
trimDir="{}/Trim".format(fastqDir)
alignDir="1.1.Align"

def getfq(wildcards):
	return samples.loc[wildcards.sampleId, ["Fq1","Fq2"]]

rule all:
	input:
		expand("1.3.Align.dedup/{sampleId}.filtered.dedup.bam", sampleId=samples.index.values.tolist()),
#		expand("1.2.Align.filtered/BaseFreq/{sampleId}.R{read}.freq.line.png", sampleId=samples.index.values.tolist(), read=[1,2])


rule trim_pe:
	input:
		getfq
	output:
		fq1="0.Fastq/Trim/{sampleId}_1.trim.fq.gz",
		fq2="0.Fastq/Trim/{sampleId}_2.trim.fq.gz"
	message:
		"Trimming... [{wildcards.sampleId}]"
	shell:
		"touch {output.fq1} {output.fq2}"

rule align_pe:
	input:
		fq1="0.Fastq/Trim/{sampleId}_1.trim.fq.gz",
		fq2="0.Fastq/Trim/{sampleId}_2.trim.fq.gz"
	output:
		"1.1.Align/{sampleId}.bam"
	message:
		"Aligning... [{wildcards.sampleId}]"
	shell:
		"touch {output}"

rule filter_align:
	input:
		"1.1.Align/{sampleId}.bam"
	output:
		"1.2.Align.filtered/{sampleId}.filtered.bam"
	message:
		"Filtering... [{wildcards.sampleId}]"
	shell:
		"touch {output}"

rule dedup_align:
	input:
		"1.2.Align.flitered/{sampleId}.filtered.bam"
	output:
		"1.3.Align.dedup/{sampleId}.filtered.dedup.bam"
	message:
		"Filtering... [{wildcards.sampleId}]"
	shell:
		"touch {output}"

#rule check_baseFreq:
#	input:
#		"1.2.Align.flitered/{sampleId}.filtered.bam"
#	output:
#		"1.2.Align.filtered/BaseFreq/{sampleId}.R1.freq.line.png",
#		"1.2.Align.filtered/BaseFreq/{sampleId}.R2.freq.line.png"
#	message:
#		"Checking baseFrequency... [{wildcards.sampleId}]"
##	shell:
#		"touch {output}"
