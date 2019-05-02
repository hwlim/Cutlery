import pandas as pd

samples = pd.read_csv("sample.tsv", sep="\t", comment="#").set_index("Id")

fastqDir="0.Fastq"
trimDir="{}/Trim".format(fastqDir)

def getfq(wildcards):
	return samples.loc[wildcards.id, ["Fq1","Fq2"]]


rule all:
	input:
		expand("0.Fastq/{id}_1.trim.fq.gz", id=samples.index.values.tolist()),
		expand("0.Fastq/{id}_2.trim.fq.gz", id=samples.index.values.tolist())

rule trim_pe:
	input:
		getfq
	output:
		fq1="{id}_1.trim.fq.gz",
		fq2="{id}_2.trim.fq.gz"
	message:
		"Trimming...{id}"
	shell:
		"touch {output.fq1} {output.fq2}"
