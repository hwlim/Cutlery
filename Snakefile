
########################
## Sample Information
import pandas as pd
samples = pd.read_csv("sample.tsv", sep="\t", comment="#")
samples_indexById = samples.set_index("Id")
#samples_indexByName = samples.set_index("Name")

########################
#fastqDir="0.Fastq"
#trimDir="{}/Trim".format(fastqDir)
#alignDir="1.1.Align"

def getfq(wildcards):
	return "0.Fastq/" + samples_indexById.loc[wildcards.sampleId, ["Fq1","Fq2"]]

#def getfq_trim(wildcards):
#	sampleId = samples_indexByName.loc[wildcards.sampleName, ["Id"]]
#	return [ "0.Fastq/Trim/" + sampleId + "_1.trim.fq.gz", "0.Fastq/Trim/" + sampleId + "_2.trim.fq.gz" ]

	
rule all:
	input:
		expand("1.2.Align.filtered/BaseFreq/{sampleId}.R{read}.freq.line.png",
			sampleId=samples["Id"].tolist(), read=[1,2]),
		expand("2.BigWig/{sampleId}.filtered.dedup.{group}.ctr.bw",
			sampleId=samples["Id"].tolist(), group=["nfr","nuc"])
#		expand("1.4.Align.split/{sampleId}.filtered.dedup.{group}.{proctype}.bed", sampleId=samples["Id"].tolist(), group=["nfr","nuc"], proctype=["con","ctr","sep"])
#		expand("1.3.Align.dedup/{sampleId}.filtered.dedup.bam", sampleId=samples["Id"].tolist()),
#		expand("1.2.Align.filtered/{sampleId}.filtered.bam", sampleId=samples.index.values.tolist()),


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
		"1.2.Align.filtered/{sampleId}.filtered.bam"
	output:
		"1.3.Align.dedup/{sampleId}.filtered.dedup.bam"
	message:
		"Deduplicating... [{wildcards.sampleId}]"
	shell:
		"touch {output}"

rule check_baseFreq:
	input:
		"1.2.Align.filtered/{sampleId}.filtered.bam"
	output:
		expand("1.2.Align.filtered/BaseFreq/{{sampleId}}.{read}.freq.line.png",
			read=["R1","R2"])
	message:
		"Checking baseFrequency... [{wildcards.sampleId}]"
	shell:
		"touch {output}"

rule split_bam:
	input:
		"1.3.Align.dedup/{sampleId}.filtered.dedup.bam"
	output:
		expand("1.4.Align.split/{{sampleId}}.filtered.dedup.{group}.{proctype}.bed",
			group=["nfr","nuc"], proctype=["con","ctr","sep"])
	message:
		"Splitting BAM file by fragment size... [{wildcards.sampleId}]"
	shell:
		"touch {output}"

rule make_bigwig:
	input:
		expand("1.4.Align.split/{{sampleId}}.filtered.dedup.{group}.ctr.bed",
			group=["nfr","nuc"])
	output:
		expand("2.BigWig/{{sampleId}}.filtered.dedup.{group}.ctr.bw",
			group=["nfr","nuc"])
	message:
		"Making bigWig files... [{wildcards.sampleId}]"
	shell:
		"touch {output}"

#rule peak_call:
#	input:
#	output:
#	message:
#	shell:
#
