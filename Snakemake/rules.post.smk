########$$$$$$$###################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
########$$$$$$$###################################
'''
Required Variables

- bamDir
- baseFreqDir
- splitDir
- fclDir
- fragLenDir
- fragAcorDir
- homerDir
- spikeinCntDir
- bigWigDir / bigWigDir1bp / 
- bigWigScaledDir / bigWigScaledDir_sub
- peak_mask
'''


rule check_baseFreq:
	input:
		bamDir + "/{sampleName}.bam"
	output:
		read1 = baseFreqDir + "/{sampleName}.R1.freq.line.png",
		read2 = baseFreqDir + "/{sampleName}.R2.freq.line.png"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		bamToBed.separate.sh -o {baseFreqDir} {input}
		checkBaseFreq.plot.sh -g {genomeFa} -o {baseFreqDir} {baseFreqDir}/{wildcards.sampleName}.filtered.R1.bed.gz
		checkBaseFreq.plot.sh -g {genomeFa} -o {baseFreqDir} {baseFreqDir}/{wildcards.sampleName}.filtered.R2.bed.gz
		"""


## Fragment Center / Length file
rule make_fcl_file:
	input:
		splitDir + "/{sampleName}.all.con.bed.gz"
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
		splitDir + "/{sampleName}.all.con.bed.gz"
		#fragDir + "/{sampleName}.frag.bed.gz"
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
		cnr.drawAutoCorFrag.r -o {fragAcorDir}/{wildcards.sampleName} {input}
		"""

rule count_spikein:
	input:
		bamDir + "/{sampleName}.bam"
#		filteredDir + "/{sampleName}.filtered.bam"
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
		bamDir + "/{sampleName}.bam"
#		get_bam_input
#		dedupDir + "/{sampleName}.dedup.bam" if doDedup else filteredDir + "/{sampleName}.filtered.bam"
	output:
		expand(splitDir + "/{{sampleName}}.{group}.{proctype}.bed.gz",
			group=["all","nfr","nuc"], proctype=["con","ctr","sep"])
	message:
		"Splitting BAM file by fragment size... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.splitBamToBed.sh -o {splitDir}/{wildcards.sampleName} -c "{chrRegexTarget}" {input}
		"""

rule make_bigwig:
	input:
		all = splitDir + "/{sampleName}.all.ctr.bed.gz",
		nfr = splitDir + "/{sampleName}.nfr.ctr.bed.gz",
		nuc = splitDir + "/{sampleName}.nuc.ctr.bed.gz"
	output:
		all = bigWigDir + "/{sampleName}.all.ctr.bw",
		nfr = bigWigDir + "/{sampleName}.nfr.ctr.bw",
		nuc = bigWigDir + "/{sampleName}.nuc.ctr.bw"
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

rule make_bigwig1bp:
	input:
		splitDir + "/{sampleName}.all.sep.bed.gz"
#		dedupDir + "/{sampleName}.dedup.bam" if doDedup else filteredDir + "/{sampleName}.filtered.bam"
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

'''
rule make_bigwig_allfrag:
	input:
		splitDir + "/{sampleName}.all.con.bed.gz"
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
'''

rule make_tagdir:
	input:
		all=splitDir + "/{sampleName}.all.ctr.bed.gz",
		nfr=splitDir + "/{sampleName}.nfr.ctr.bed.gz",
		nuc=splitDir + "/{sampleName}.nuc.ctr.bed.gz"
	output:
		all=directory(homerDir + "/{sampleName}/TSV.all"),
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
		cnr.makeHomerDir.sh -o {output.all} -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.all}
		"""



#def get_ctrl(wildcards):
#	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
#	if ctrlName.tolist()[0].upper() == "NULL":
#		return "NULL"
#	else:
#		return homerDir + "/" + ctrlName + "/TSV"

## Returns peak calling input tagDir(s): ctrl (optional) & target
def get_peakcall_input(sampleName, fragment):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ homerDir + "/" + sampleName + "/TSV." + fragment ]
	else:
		return [ homerDir + "/" + ctrlName + "/TSV." + fragment, homerDir + "/" + sampleName + "/TSV." + fragment ]

'''
def get_peakcall_factor_input(wildcards):
	# return ordered [ctrl , target] list. if no ctrl, simply [target].
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ homerDir + "/" + wildcards.sampleName + "/TSV.nfr" ]
	else:
		return [ homerDir + "/" + ctrlName + "/TSV.nfr", homerDir + "/" + wildcards.sampleName + "/TSV.nfr" ]
'''

rule call_peaks_factor:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName,"nfr")
#		get_peakcall_factor_input
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

'''
def get_peakcall_factor_input_allfrag(wildcards):
	# return ordered [ctrl , target] list. if no ctrl, simply [target].
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ homerDir + "/" + wildcards.sampleName + "/TSV.all" ]
	else:
		return [ homerDir + "/" + ctrlName + "/TSV.all", homerDir + "/" + wildcards.sampleName + "/TSV.all" ]
'''

rule call_peaks_factor_allfrag:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName,"all")
#		get_peakcall_factor_input_allfrag
	output:
		homerDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.homer.exBL.1rpm.bed"
	params:
		mask = peak_mask,
		peakDir = homerDir + "/{sampleName}/HomerPeak.factor.allFrag",
		optStr = lambda wildcards, input: "-i" if len(input)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load CnR/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {params.mask} -s \"-fragLength 100\" {params.optStr} {input}
		"""

'''
def get_peakcall_histone_input(wildcards):
	# return ordered [ctrl , target] list. if no ctrl, simply [target].
	ctrlName = samples.Ctrl[samples.Name == wildcards.sampleName]
	ctrlName = ctrlName.tolist()[0]
	if ctrlName.upper() == "NULL":
		return [ homerDir + "/" + wildcards.sampleName + "/TSV.nuc" ]
	else:
		return [ homerDir + "/" + ctrlName + "/TSV.nuc", homerDir + "/" + wildcards.sampleName + "/TSV.nuc" ]
'''

rule call_peaks_histone:
	input:
		lambda wildcards: get_peakcall_input(wildcards.sampleName,"nuc")
#		get_peakcall_histone_input
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

#### Note: Rules below simply copied from ChIP-seq rules. Needs revision for CUT&Run
'''

rule make_bigwig_scaled_divide:
	input:
		get_bigwig_scaled_input
	output:
		bigWigScaledDir_div + "/{sampleName}.scaled.divInput.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	#params:
	#	memory = "%dG" % (  cluster["make_bigwig_subtract"]["memory"]/1000 - 1 )
	shell:
		"""
		module load CnR/1.0
		bigWigDivide.sh -g {chrom_size} -m 5G -s log -a 1 -o {output} {input}
		"""




##### Average bigWig files
def get_bigwig_rep_sub(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigDir_sub + "/" + x + ".subInput.bw", repL)

rule make_bigwig_sub_avg:
	input:
		get_bigwig_rep_sub
	output:
		bigWigDir_sub_avg + "/{groupName}.subInput.avg.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""

def get_bigwig_rep_scaled_sub(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigScaledDir_sub + "/" + x + ".subInput.bw", repL)

rule make_bigwig_scaled_sub_avg:
	input:
		get_bigwig_rep_scaled_sub
	output:
		bigWigScaledDir_sub_avg + "/{groupName}.scaled.subInput.avg.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""

def get_bigwig_rep_scaled_div(wildcards):
	repL = samples.Name[samples.Group == wildcards.groupName].tolist()
	return map(lambda x: bigWigScaledDir_div + "/" + x + ".divInput.bw", repL)

rule make_bigwig_scaled_div_avg:
	input:
		get_bigwig_rep_scaled_div
	output:
		bigWigScaledDir_div_avg + "/{groupName}.scaled.divInput.avg.bw"
	message:
		"Making average bigWig files... [{wildcards.groupName}]"
	params:
		memory = "5G"
	shell:
		"""
		module load CnR/1.0
		makeBigWigAverage.sh -g {chrom_size} -m {params.memory} -o {output} {input}
		"""
'''

