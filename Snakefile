SRCFILES = ["1", "2", "3"]

rule all:
	input:
		"sortpasted.txt"

rule sort:
	input:
		"src{sample}.txt"
	output:
		"src{sample}.sorted.txt"
	shell:
		"sort -k1,1n {input} > {output}"

rule paste:
	input:
		expand("src{sample}.sorted.txt", sample=SRCFILES)
	output:
		"sortpasted.txt"
	shell:
		"paste {input} > {output}"

