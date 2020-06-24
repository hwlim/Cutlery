################################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
################################################



#################################
## Default values for undefined variables in Snakefile

## bamDir is used to designated bam file directory for downstream analysis
## Originally motivated for reuse or rule.post.smk in replicate pooling in Snakemake.pool
## If bamDir is not defined, i.e. for simply processing individual replicates not pooling
## dedupDir or filteredDir is selected 
if "bamDir" not in locals():
	if doDedup:
		bamDir = dedupDir
	else:
		bamDir = filteredDir



#################################
## Rule Start


## Nucleotide base frequence around 5'-ends of read 1 & 2
## NEED REVISION for OUTPUT names & directories
rule check_baseFreq:
	input:
		bamDir + "/{sampleName}.bam"
	output:
		sampleDir + "/{sampleName}/QC/base_freq.R1.freq.png",
		sampleDir + "/{sampleName}/QC/base_freq.R2.freq.png"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		bamToBed.separate.sh -o {sampleDir}/{wildcards.sampleName}/tmp {input}
		checkBaseFreq.plot.sh -g {genomeFa} -n {wildcards.sampleName} -o {sampleDir}/{wildcards.sampleName}/QC/base_freq.R1 {sampleDir}/{wildcards.sampleName}/tmp.R1.bed.gz
		checkBaseFreq.plot.sh -g {genomeFa} -n {wildcards.sampleName} -o {sampleDir}/{wildcards.sampleName}/QC/base_freq.R2 {sampleDir}/{wildcards.sampleName}/tmp.R2.bed.gz
		rm {sampleDir}/{wildcards.sampleName}/tmp.R1.bed.gz
		rm {sampleDir}/{wildcards.sampleName}/tmp.R2.bed.gz
		"""


## BAM to fragment bed files: all / nfr / nuc
rule split_bam:
	input:
		bamDir + "/{sampleName}.bam"
	output:
		expand(sampleDir + "/{{sampleName}}/Fragments/frag.{group}.{proctype}.bed.gz",
			group=["all","nfr","nuc"], proctype=["con","ctr","sep"])
	message:
		"Splitting BAM file by fragment size... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.splitBamToBed.sh -o {sampleDir}/{wildcards.sampleName}/Fragments/frag -c "{chrRegexTarget}" {input}
		"""


## FCL (Fragment Center / Length) file
rule make_fcl_file:
	input:
		sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz"
	output:
		sampleDir + "/{sampleName}/fcl.bed.gz"
#	params:
#		memory = "%dG" % ( cluster["make_fcl_file"]["memory"]/1000 - 1 )
	message:
		"Making FCL bed files... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		fragmentToFCL.sh -o {output} -m 5G {input}
		"""

## Fragment length distribution
rule get_fragLenHist:
	input:
		sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/fragLen.dist.txt",
		sampleDir + "/{sampleName}/QC/fragLen.dist.png"
	message:
		"Checking fragment length... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		ngs.fragLenHist.r -o {sampleDir}/{wildcards.sampleName}/QC/fragLen -n {wildcards.sampleName} {input}
		"""

## Autocorrelation plot as Q/C
rule get_frag_autocor:
	input:
		sampleDir + "/{sampleName}/fcl.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/enrich.acor.txt",
		sampleDir + "/{sampleName}/QC/enrich.acor.png"
	message:
		"Checking fragment auto-correlation... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.drawAutoCorFrag.r -o {sampleDir}/{wildcards.sampleName}/QC/enrich {input}
		"""

rule count_spikein:
	input:
		bamDir + "/{sampleName}.bam"
	output:
		sampleDir + "/{sampleName}/QC/spikeCnt.txt"
	message:
		"Counting spikein tags... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		countSpikein.sh -p {spikePrefix} {input} > {output}
		"""

rule make_spikeintable:
	input:
		expand(sampleDir + "/{sampleName}/QC/spikeCnt.txt", sampleName=samples.Name.tolist())
	output:
		qcDir + "/spikein.txt"
	message:
		"Making spikein table..."
	shell:
		"""
		module load CnR/1.0
		makeSpikeCntTable.r -o {spikeinCntDir}/spikein {input}
		"""

rule make_bigwig:
	input:
		all = sampleDir + "/{sampleName}/Fragments/frag.all.ctr.bed.gz",
		nfr = sampleDir + "/{sampleName}/Fragments/frag.nfr.ctr.bed.gz",
		nuc = sampleDir + "/{sampleName}/Fragments/frag.nuc.ctr.bed.gz"
	output:
		all = sampleDir + "/{sampleName}/igv.all.ctr.bw",
		nfr = sampleDir + "/{sampleName}/igv.nfr.ctr.bw",
		nuc = sampleDir + "/{sampleName}/igv.nuc.ctr.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -o {output.all} {input.all}
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -o {output.nfr} {input.nfr}
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -o {output.nuc} {input.nuc}
		"""


## 1bp-resolution bigwig files
## NOTE: considering protrusion handling, it's better to use fragment file than seprate reads
## 		Thus needs an update in command line
rule make_bigwig1bp:
	input:
		sampleDir + "/{sampleName}/Fragments/frag.all.sep.bed.gz"
	output:
		sampleDir + "/{sampleName}/igv.1bp.plus.bw",
		sampleDir + "/{sampleName}/igv.1bp.minus.bw"
	message:
		"Making 1bp-resolution bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		ngs.alignToBigWig.sh -o {sampleDir}/{wildcards.sampleName}/igv.1bp -g {chrom_size} -l 1 -m 5G -c "{chrRegexTarget}" {input}
		"""

rule make_bigwig_allfrag:
	input:
		sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz"
	output:
		sampleDir + "/{sampleName}/igv.allFrag.bw"
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
		all = sampleDir + "/{sampleName}/Fragments/frag.all.ctr.bed.gz",
		nfr = sampleDir + "/{sampleName}/Fragments/frag.nfr.ctr.bed.gz",
		nuc = sampleDir + "/{sampleName}/Fragments/frag.nuc.ctr.bed.gz"
	output:
		all=directory(sampleDir + "/{sampleName}/TSV.all"),
		nfr=directory(sampleDir + "/{sampleName}/TSV.nfr"),
		nuc=directory(sampleDir + "/{sampleName}/TSV.nuc")
#	params:
#		name = "{sampleName}"
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.makeHomerDir.sh -o {output.nfr} -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.nfr}
		cnr.makeHomerDir.sh -o {output.nuc} -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.nuc}
		cnr.makeHomerDir.sh -o {output.all} -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.all}
		"""



## Returns peak calling input tagDir(s): ctrl (optional) & target
def get_peakcall_input(sampleName, fragment):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ sampleDir + "/" + sampleName + "/TSV." + fragment ]
	else:
		return [ sampleDir + "/" + ctrlName + "/TSV." + fragment, sampleDir + "/" + sampleName + "/TSV." + fragment ]


rule call_peaks_factor:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName,"nfr")
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	params:
		mask = peak_mask,
		peakDir = sampleDir + "/{sampleName}/HomerPeak.factor",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""


rule call_peaks_factor_allfrag:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName,"all")
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed"
	params:
		mask = peak_mask,
		peakDir = sampleDir + "/{sampleName}/HomerPeak.factor.allFrag",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""


rule call_peaks_histone:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName,"nuc")
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed"
	params:
		mask = peak_mask,
		peakDir = sampleDir + "/{sampleName}/HomerPeak.histone",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.peakCallHistone.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""


rule call_peaks_histone_allfrag:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName,"all")
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone.allFrag/peak.exBL.bed"
	params:
		mask = peak_mask,
		peakDir = sampleDir + "/{sampleName}/HomerPeak.histone.allFrag",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.peakCallHistone.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""

rule run_homermotif:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed.all.noBG/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load Motif/1.0
		runHomerMotif.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse -o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor {input}
		"""


rule draw_peak_heatmap_factor:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/heatmap.exBL.1rpm.png"
	params:
		workDir = sampleDir + "/{sampleName}"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		drawBigWigHeatmap.r -t {wildcards.sampleName} -m 0,0.5,2,0.5 -w 2000 -c NFR,NUC -s 3,6 \
			-o {params.workDir}/HomerPeak.factor/heatmap.exBL.1rpm \
			{input} {params.workDir}/igv.nfr.ctr.bw {params.workDir}/igv.nuc.ctr.bw
		"""

#####################################################
## Scaled BigWig by Spike-in using raw read counts (not RPM)
#def get_scalefactor(wildcards):
#	# return ordered [ctrl , target] list.
#	spikeinTable = pd.read_csv(spikeinCntDir + "/spikein.txt", sep="\t", comment="#", na_filter=False)
#	if not spikeinTable.Sample.is_unique:
#		print( "Error: Sample column spikein.txt is not unique")
#		sys.exit()
#	return [ spikeinTable.ScaleFactor[spikeinTable.Sample == wildcards.sampleName].tolist()[0] ]

################################################################
## *** Note *****
## get_scalefactor function is deterministic, i.e Snakemake run this function to generate command lines to submit
## Therefore, when spikein.txt does not exist, this rule will invoke error.
## This rule must be revised the command take spikein.txt as an input directly not using the get_scalefactor function
#  

'''
## NOTE: this rule is deprecate. planning to use RPSM below
## Raw read count scale + spike-in scaled
rule make_bigwig_scaled:
	input:
		all = splitDir + "/{sampleName}.all.ctr.bed.gz",
		nfr = splitDir + "/{sampleName}.nfr.ctr.bed.gz",
		nuc = splitDir + "/{sampleName}.nuc.ctr.bed.gz",
		spikeinCnt = spikeinCntDir + "/spikein.txt"
	output:
		all = bigWigScaledDir + "/{sampleName}.all.ctr.bw",
		nfr = bigWigScaledDir + "/{sampleName}.nfr.ctr.bw",
		nuc = bigWigScaledDir + "/{sampleName}.nuc.ctr.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "%dG" % ( cluster["make_bigwig"]["memory"]/1000 - 1 ),
#		scaleFactor = get_scalefactor
	shell:
		"""
		module load CnR/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | cut -f 6`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.all} {input.all}
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.nfr} {input.nfr}
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.nuc} {input.nuc}
		"""


## NOTE: this rule is deprecate. planning to use RPSM below. needs revision for RPSM
def get_bigwig_scaled_input(sampleName, fragment):
	# return ordered [ctrl , target] list.
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	return [ bigWigScaledDir + "/" + sampleName + "." + fragment + ".ctr.bw",
			bigWigScaledDir + "/" + ctrlName + "." + fragment + ".ctr.bw"]

rule make_bigwig_scaled_subtract:
	input:
		all=lambda wildcards: get_bigwig_scaled_input(wildcards.sampleName,"all"),
		nfr=lambda wildcards: get_bigwig_scaled_input(wildcards.sampleName,"nfr"),
		nuc=lambda wildcards: get_bigwig_scaled_input(wildcards.sampleName,"nuc")
	output:
		all=bigWigScaledDir_sub + "/{sampleName}.all.scaled.subInput.bw",
		nfr=bigWigScaledDir_sub + "/{sampleName}.nfr.scaled.subInput.bw",
		nuc=bigWigScaledDir_sub + "/{sampleName}.nuc.scaled.subInput.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output.all} {input.all}
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output.nfr} {input.nfr}
		bigWigSubtract.sh -g {chrom_size} -m 5G -t -1000 {output.nuc} {input.nuc}
		"""


rule make_bigwig_allfrag_rpsm:
	input:
		bed = sampleDir + "/{sampleName}/Fragments/frag.all.ctr.bed.gz",
		sampleDir + "/{sampleName}/QC/spikeCnt.txt"
	output:
		sampleDir + "/{sampleName}/igv.rpsm.allFrag.bw"
	message:
		"Making spike-in scaled allFrag bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | gawk '{{ printf "%f", 100000/$3 }}'`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		cnr.bedToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""

'''