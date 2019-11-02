########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
########$$$$$$$###################################

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
		filteredDir + "/{sampleName}.filtered.bam"
	message:
		"Filtering... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.filterBam.sh  -o {output} -c "{chrRegexAll}" {input}
		"""

rule dedup_align:
	input:
		filteredDir + "/{sampleName}.filtered.bam"
	output:
		dedupDir + "/{sampleName}.dedup.bam"
	message:
		"Deduplicating... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["dedup_align"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		cnr.dedupBam.sh -m {params.memory} -o {output} -r {input}
		"""

rule check_baseFreq:
	input:
		filteredDir + "/{sampleName}.filtered.bam"
	output:
		read1 = baseFreqDir + "/{sampleName}.filtered.R1.freq.line.png",
		read2 = baseFreqDir + "/{sampleName}.filtered.R2.freq.line.png"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		bamToBed.separate.sh -o {baseFreqDir} {input}
		checkBaseFreq.plot.sh -g {genomeFa} -o {baseFreqDir} {baseFreqDir}/{wildcards.sampleName}.filtered.R1.bed.gz
		checkBaseFreq.plot.sh -g {genomeFa} -o {baseFreqDir} {baseFreqDir}/{wildcards.sampleName}.filtered.R2.bed.gz
		"""

## fragment bed file for V-plot analysis ** NOT FOR PEAK-CALLING
rule make_fragment:
	input:
		dedupDir + "/{sampleName}.dedup.bam" if doDedup else filteredDir + "/{sampleName}.filtered.bam"
	output:
		fragDir + "/{sampleName}.frag.bed.gz"
#	params:
#		memory = "%dG" % ( cluster["make_fragment"]["memory"]/1000 - 1 )
	message:
		"Making fragment bed files... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		bamToFragment.sh -o {output} -l -1 -s -m 5G {input}
		"""

rule make_fcl_file:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		fclDir + "/{sampleName}.fcl.bed.gz"
#	params:
#		memory = "%dG" % ( cluster["make_fcl_file"]["memory"]/1000 - 1 )
	message:
		"Making FCL bed files... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		fragmentToFCL.sh -o {output} -m 5G {input}
		"""

rule get_fragLenHist:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		fragLenDir + "/{sampleName}.dist.txt",
		fragLenDir + "/{sampleName}.dist.png"
	message:
		"Checking fragment length... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		ngs.fragLenHist.r -o {fragLenDir}/{wildcards.sampleName} {input}
		"""

rule get_frag_autocor:
	input:
		fclDir + "/{sampleName}.fcl.bed.gz"
	output:
		fragAcorDir + "/{sampleName}.acor.txt",
		fragAcorDir + "/{sampleName}.acor.png"
	message:
		"Checking fragment auto-correlation... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.drawAutoCorFrag.r -o {fragAcorDir}/{wildcards.sampleName} -d 1000 {input}
		"""

rule count_spikein:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		spikeinCntDir + "/{sampleName}.spikeCnt.txt"
	message:
		"Counting spikein tags... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		countSpikein.sh -p {spikePrefix} {input} > {output}
		"""

rule make_spikeintable:
	input:
		expand(spikeinCntDir + "/{sampleName}.spikeCnt.txt", sampleName=samples.Name.tolist())
	output:
		spikeinCntDir + "/spikein.txt"
	message:
		"Making spikein table..."
	shell:
		"""
		module load CnR/1.0
		makeSpikeCntTable.r -o {spikeinCntDir}/spikein {input}
		"""


rule split_bam:
	input:
		dedupDir + "/{sampleName}.dedup.bam" if doDedup else filteredDir + "/{sampleName}.filtered.bam"
	output:
		expand(splitDir + "/{{sampleName}}.{group}.{proctype}.bed.gz",
			group=["nfr","nuc"], proctype=["con","ctr","sep"])
	message:
		"Splitting BAM file by fragment size... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.splitBamToBed.sh -o {splitDir}/{wildcards.sampleName} -c "{chrRegexTarget}" {input}
		"""

rule make_bigwig:
	input:
		nfr = splitDir + "/{sampleName}.nfr.ctr.bed.gz",
		nuc = splitDir + "/{sampleName}.nuc.ctr.bed.gz"
	output:
		nfr = bigWigDir + "/{sampleName}.nfr.ctr.bw",
		nuc = bigWigDir + "/{sampleName}.nuc.ctr.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -o {output.nfr} {input.nfr}
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -o {output.nuc} {input.nuc}
		"""

rule make_bigwig1bp:
	input:
		dedupDir + "/{sampleName}.dedup.bam" if doDedup else filteredDir + "/{sampleName}.filtered.bam"
	output:
		bigWigDir1bp + "/{sampleName}.plus.bw",
		bigWigDir1bp + "/{sampleName}.minus.bw"
	message:
		"Making 1bp-resolution bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		ngs.alignToBigWig.sh -o {bigWigDir1bp}/{wildcards.sampleName} -g {chrom_size} -l 1 -m 5G -c "{chrRegexTarget}" {input}

		"""

rule make_bigwig_allfrag:
	input:
		fragDir + "/{sampleName}.frag.bed.gz"
	output:
		bigWigDirAllFrag + "/{sampleName}.allFrag.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -o {output} {input}
		"""


rule make_tagdir:
	input:
		nfr=splitDir + "/{sampleName}.nfr.ctr.bed.gz",
		nuc=splitDir + "/{sampleName}.nuc.ctr.bed.gz"
	output:
		nfr=directory(homerDir + "/{sampleName}/TSV.nfr"),
		nuc=directory(homerDir + "/{sampleName}/TSV.nuc")
#	params:
#		name = "{sampleName}"
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.makeHomerDir.sh -o {output.nfr} -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.nfr}
		cnr.makeHomerDir.sh -o {output.nuc} -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.nuc}
		"""



#def get_ctrl(wildcards):
#	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
#	if ctrlName.tolist()[0].upper() == "NULL":
#		return "NULL"
#	else:
#		return homerDir + "/" + ctrlName + "/TSV"

def get_peakcall_factor_input(wildcards):
	# return ordered [ctrl , target] list. if no ctrl, simply [target].
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ homerDir + "/" + wildcards.sampleName + "/TSV.nfr" ]
	else:
		return [ homerDir + "/" + ctrlName + "/TSV.nfr", homerDir + "/" + wildcards.sampleName + "/TSV.nfr" ]

rule call_peaks_factor:
	input:
		get_peakcall_factor_input
	output:
		homerDir + "/{sampleName}/HomerPeak.factor/peak.homer.exBL.1rpm.bed"
	params:
		mask = peak_mask,
		peakDir = homerDir + "/{sampleName}/HomerPeak.factor",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""

def get_peakcall_histone_input(wildcards):
	# return ordered [ctrl , target] list. if no ctrl, simply [target].
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ homerDir + "/" + wildcards.sampleName + "/TSV.nuc" ]
	else:
		return [ homerDir + "/" + ctrlName + "/TSV.nuc", homerDir + "/" + wildcards.sampleName + "/TSV.nuc" ]

rule call_peaks_histone:
	input:
		get_peakcall_histone_input
	output:
		homerDir + "/{sampleName}/HomerPeak.histone/peak.homer.exBL.bed"
	params:
		mask = peak_mask,
		peakDir = homerDir + "/{sampleName}/HomerPeak.histone",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.peakCallHistone.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""


#####################################################
## Rules for data pooling by "Group" column

def get_bigwig_rep_nfr(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigDir + "/" + x + ".nfr.ctr.bw", repL)

def get_bigwig_rep_nuc(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigDir + "/" + x + ".nuc.ctr.bw", repL)

rule make_bigwig_avg:
	input:
		nfr = get_bigwig_rep_nfr,
		nuc = get_bigwig_rep_nuc
	output:
		nfr = bigWigDir_avg + "/{groupName}.nfr.ctr.bw",
		nuc = bigWigDir_avg + "/{groupName}.nuc.ctr.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		makeBigWigAverage.sh -g {chrom_size} -m 5G -o {output.nfr} {input.nfr}
		makeBigWigAverage.sh -g {chrom_size} -m 5G -o {output.nuc} {input.nuc}
		"""

