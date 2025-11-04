################################################
### Snakemake rules for CUT&RUN Data Processing
### - except for "rule all"
### 
###		Written by Hee Woong Lim
################################################


################################################
## TODO:
##	1. required file validataion
##	2, required variables validation
##	3. MEME motif search & footprinting


#################################
## Default values for undefined variables in Snakefile
import sys


if "meme_db" not in locals():
	meme_db = os.environ["LIMLAB_BASE"] + "/Motif/MEME_DB/Merged_By_Lim.meme"

if "numHighestPeaks" not in locals():
	numHighestPeaks = 5

## 2kb bed file of all promoters for ATAC-seq QC
if "bed_promoter" not in locals():
	bed_promoter = "NULL"

if 'species_macs' not in locals():
	#print("Warning: species_macs is not defined; using hs (default)", file=sys.stderr)
	species_macs="hs"

if 'do_csem' not in locals():
	#print("Warning: do_csem is not defined; setting False", file=sys.stderr)
	do_csem=False

#################################
## Rule Start

# ## bamDir is used to designated bam file directory for downstream analysis
# ## Originally motivated for reuse or rule.post.smk in replicate pooling in Snakemake.pool
# ## If bamDir is not defined, i.e. for simply processing individual replicates not pooling
# ## dedupDir or alignDir is selected 
# if "bamDir" not in locals():
# 	if doDedup:
# 		bamDir = dedupDir
# 	else:
# 		bamDir = alignDir

## 20240710: Need update to handle CSEM results
rule dedup_align:
	input:
		alignDir + "/{sampleName}/align.bam"
	output:
		bam = dedupDir + "/{sampleName}/align.bam",
		bai = dedupDir + "/{sampleName}/align.bam.bai"
	message:
		"Deduplicating... [{wildcards.sampleName}]"
	params:
		memory = "%dG" % ( cluster["dedup_align"]["memory"]/1000 - 1 )
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.dedupBam.sh -m {params.memory} -o {output.bam} -r {input}
		samtools index {output.bam}
		"""


## find bam file folder for single-end style homer tag directory creation
def get_bam_for_downstream(sampleName):
	# return ordered [ctrl , target] list.
	if doDedup:
		bam = dedupDir + "/" + sampleName + "/align.bam"
	else:
		if do_csem:
			bam = alignDir + "/" + sampleName + "/CSEM/align.uniq.bam"
		else:
			bam = alignDir + "/" + sampleName + "/align.bam"

	return bam


## Convert BAM to fragment bed file
## - Chromosome filtering
## - 0x2  : Concordant pairs only
## - 0x400: Removes duplicates
rule make_fragment:
	input:
		bam = lambda wildcards: get_bam_for_downstream(wildcards.sampleName),
		bai = lambda wildcards: get_bam_for_downstream(wildcards.sampleName) + ".bai"
	output:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	message:
		"Making fragment bed files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		mkdir -p {sampleDir}/{wildcards.sampleName}
		ngs.bamToFragment.py -c "{chrRegexAll}" -f 0x2 -F 0x400 {input.bam} | sort -S 4G -k1,1 -k2,2n -k3,3n | gzip > {output}
		"""


## Count the number of unique fragents out of total
## - Number & % of unique fragments
rule count_uniq_fragment:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/fragment.uniq_cnt.txt"
	message:
		"Counting unique fragments... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.countUniqFrag.sh -n {wildcards.sampleName} -o {output} {input}
		"""


## Generate unique fragment count table by combining information from all sample
## And generate diagnostic plots
rule make_uniqcnt_table:
	input:
		expand(sampleDir + "/{sampleName}/QC/fragment.uniq_cnt.txt", sampleName=samples.Name.tolist())
	output:
		expand(qcDir + "/uniqFragCnt.{ext}", ext=[ "txt", "pdf", "png" ])
	message:
		"Making unique fragment count table..."
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.makeUniqCntTable.r -o {qcDir}/uniqFragCnt {input}
		"""

## Nucleotide base frequence around 5'-ends of read 1 & 2
## NEED REVISION for OUTPUT names & directories
rule check_baseFreq:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		genomeFa = genomeFa
	output:
		sampleDir + "/{sampleName}/QC/base_freq.png",
		sampleDir + "/{sampleName}/QC/base_freq.html"

	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		checkBaseFreq.plot.sh -o {sampleDir}/{wildcards.sampleName}/QC/base_freq \
			-n {wildcards.sampleName} -g {input.genomeFa} -c "{chrRegexTarget}" -m both -l 20 -f -i -v {input.frag}
		"""


rule check_baseFreq_chrM:
	input:
		bam = lambda wildcards: get_bam_for_downstream(wildcards.sampleName)
	output:
		expand(sampleDir + "/{{sampleName}}/QC/base_freq_chrM.{ext}", ext=["png","html"])
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		checkBaseFreq.plot.sh -o {sampleDir}/{wildcards.sampleName}/QC/base_freq_chrM \
			-n {wildcards.sampleName} -g {genomeFa} -c "^chrM" -m 5 -l 20 -i -v {input.bam}
		"""


## BAM to fragment bed files: all / nfr / nuc
rule split_bam:
	input:
		bam = lambda wildcards: get_bam_for_downstream(wildcards.sampleName)
	output:
		expand(sampleDir + "/{{sampleName}}/Fragments/frag.{group}.{proctype}.bed.gz",
			group=["all","nfr","nuc"], proctype=["con","ctr","sep"])
	message:
		"Splitting BAM file by fragment size... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.splitBamToBed.sh -o {sampleDir}/{wildcards.sampleName}/Fragments/frag -c "{chrRegexTarget}" {input.bam}
		"""

## split fragment file by length 
## NFR / NUC
rule split_fragment_ctr:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/fragment.ctr.nfr.bed.gz",
		sampleDir + "/{sampleName}/fragment.ctr.nuc.bed.gz"
	message:
		"Splitting fragment file by fragment size... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.splitFragFile.sh -o {sampleDir}/{wildcards.sampleName}/fragment.ctr -c "{chrRegexTarget}" {input}
		"""

## FCL (Fragment Center / Length) file
rule make_fcl_file:
	input:
		#sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz"
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/fcl.bed.gz"
#	params:
#		memory = "%dG" % ( cluster["make_fcl_file"]["memory"]/1000 - 1 )
	message:
		"Making FCL bed files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		fragmentToFCL.sh -o {output} -m 5G {input}
		"""

## Fragment length distribution
## - distribution data file
## - distribution plot
rule get_fragLenHist:
	input:
		#sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz"
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/fragLen.dist.txt",
		sampleDir + "/{sampleName}/QC/fragLen.dist.png"
	message:
		"Checking fragment length... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragLenHist.r -o {sampleDir}/{wildcards.sampleName}/QC/fragLen.dist -n {wildcards.sampleName} {input}
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
		module purge
		module load Cutlery/1.0
		cnr.drawAutoCorFrag.r -o {sampleDir}/{wildcards.sampleName}/QC/enrich {input}
		"""

rule count_spikein:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/spikeCnt.txt"
	message:
		"Counting spikein tags... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.countSpikein.sh -p {spikePrefix} -n {wildcards.sampleName} {input} > {output}
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
		module purge
		module load Cutlery/1.0
		ngs.makeSpikeCntTable.r -o {qcDir}/spikein {input}
		"""

## Make raw bedgraph files from fragment.bed for SEACR input
rule make_bedgraph_frag:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		chrom = chrom_size
	output:
		all=sampleDir + "/{sampleName}/igv.raw.bedGraph.gz"
	message:
		"Making raw bedgraph files from fragment.bed... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragToBigWig.sh -b -g {input.chrom} -c "{chrRegexTarget}" -m 10G -s 1 -o {output.all} {input.frag}
		"""

## Returns peak calling input tagDir(s): ctrl (optional) & target
def get_seacr_input(sampleName):
	ctrlName = get_ctrl_name(sampleName)
	if ctrlName.upper() == "NULL":
		return [ sampleDir + "/" + sampleName + "/igv.raw.bedGraph.gz" ]
	else:
		return [ sampleDir + "/" + sampleName + "/igv.raw.bedGraph.gz", sampleDir + "/" + ctrlName + "/igv.raw.bedGraph.gz" ]

def get_seacr_param(sampleName):
	ctrlName = get_ctrl_name(sampleName)
	if ctrlName.upper() == "NULL":			
		return [ sampleDir + "/" + sampleName + "/igv.raw.bedGraph.gz", "0.01 non" ]
	else:
		return [ sampleDir + "/" + sampleName + "/igv.raw.bedGraph.gz", sampleDir + "/" + ctrlName + "/igv.raw.bedGraph.gz", "norm" ]

## Run SEACR
rule run_seacr:
	input:
		lambda wildcards: get_seacr_input(wildcards.sampleName)
	output:
		expand(sampleDir + "/{{sampleName}}/SEACR/peak.{type}.{ext}", type=["stringent", "relaxed"], ext=["exBL.bed", "bed"])
	message:
		"Running SEACR... [{wildcards.sampleName}]"
	params:
		input = lambda wildcards: get_seacr_param(wildcards.sampleName)
	shell:
		"""
		module purge
		module load Cutlery
		SEACR_1.3.sh {params.input} stringent {sampleDir}/{wildcards.sampleName}/SEACR/peak
		SEACR_1.3.sh {params.input} relaxed {sampleDir}/{wildcards.sampleName}/SEACR/peak

		## Blacklist filtering if given
		if [ {peak_mask} != "NULL" ];then
			intersectBed -a {sampleDir}/{wildcards.sampleName}/SEACR/peak.stringent.bed -b {peak_mask} -v > {sampleDir}/{wildcards.sampleName}/SEACR/peak.stringent.exBL.bed
			intersectBed -a {sampleDir}/{wildcards.sampleName}/SEACR/peak.relaxed.bed -b {peak_mask} -v > {sampleDir}/{wildcards.sampleName}/SEACR/peak.relaxed.exBL.bed
		else
			cat {sampleDir}/{wildcards.sampleName}/SEACR/peak.stringent.bed > {sampleDir}/{wildcards.sampleName}/SEACR/peak.stringent.exBL.bed
			cat {sampleDir}/{wildcards.sampleName}/SEACR/peak.relaxed.bed > {sampleDir}/{wildcards.sampleName}/SEACR/peak.relaxed.exBL.bed
		fi
		"""


## BigWig files from 100bp re-sized fragments
rule make_bigwig_ctr:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		chrom = chrom_size
	output:
		all = sampleDir + "/{sampleName}/igv.all.ctr.bw",
		nfr = sampleDir + "/{sampleName}/igv.nfr.ctr.bw",
		nuc = sampleDir + "/{sampleName}/igv.nuc.ctr.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -r 100 -m 5G -o {output.all} {input.frag}
		ngs.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 0 -L 119 -r 100 -m 5G -o {output.nfr} {input.frag}
		ngs.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 151 -L 1000000 -r 100 -m 5G -o {output.nuc} {input.frag}
		"""

## BigWig files from original sized fragments
rule make_bigwig_frag_all:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		chrom = chrom_size
	output:
		sampleDir + "/{sampleName}/igv.all.con.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -m 5G -o {output} {input.frag}
		"""

rule make_bigwig_frag_nfr:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		chrom = chrom_size
	output:
		sampleDir + "/{sampleName}/igv.nfr.con.bw",
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 0 -L 119 -m 5G -o {output} {input.frag}
		"""

rule make_bigwig_frag_nuc:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		chrom = chrom_size
	output:
		sampleDir + "/{sampleName}/igv.nuc.con.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 151 -L 1000000 -m 5G -o {output} {input.frag}
		"""

## BigWig files in 1bp-resolution using 5'-end 1bp
rule make_bigwig1bp:
	input:
		#frag = sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz",
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		chrom = chrom_size
		#sampleDir + "/{sampleName}/Fragments/frag.all.sep.bed.gz"
	output:
		sampleDir + "/{sampleName}/igv.1bp.plus.bw",
		sampleDir + "/{sampleName}/igv.1bp.minus.bw"
	message:
		"Making 1bp-resolution bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragToBigWigStranded1bp.sh -o {sampleDir}/{wildcards.sampleName}/igv.1bp -g {input.chrom} -c "{chrRegexTarget}" -s 0 -m 5G {input.frag}
		"""
#		ngs.alignToBigWig.sh -o {sampleDir}/{wildcards.sampleName}/igv.1bp -g {chrom_size} -l 1 -m 5G -c "{chrRegexTarget}" {input}


## 5'-end 1bp BigWig file in raw read count scale & nonnegative values (positive values both for plus/minus track)
rule make_bigwig1bp_raw_abs:
	input:
		#frag = sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz",
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		chrom = chrom_size
		#sampleDir + "/{sampleName}/Fragments/frag.all.sep.bed.gz"
	output:
		sampleDir + "/{sampleName}/igv.1bp.raw.abs.plus.bw",
		sampleDir + "/{sampleName}/igv.1bp.raw.abs.minus.bw"
	message:
		"Making 1bp-resolution raw/abs bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.fragToBigWigStranded1bp.sh -o {sampleDir}/{wildcards.sampleName}/igv.1bp.raw.abs -g {input.chrom} -c "{chrRegexTarget}" -s 1 -m 5G -n {input.frag}
		"""

## RNA-seq style bigwig file generation ## to correctly visualize spliced fragments
## Incorporated to handle Collins ATAC-seq fly SPS data (Brian lab)
rule make_bigwig_splice:
	input:
		bam = lambda wildcards: get_bam_for_downstream(wildcards.sampleName),
		bai = lambda wildcards: get_bam_for_downstream(wildcards.sampleName) + ".bai",
		chrom = chrom_size
	output:
		sampleDir + "/{sampleName}/igv.all.splice.bw"
	#params:
	#	paired = "-p" if is_paired("{wildcards.sampleName}") else ""
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load RNAseq/1.0
		rnaSeq.bamToBigWig.sh -g {input.chrom} -m unstranded -p -c "{chrRegexTarget}" -v -o {sampleDir}/{wildcards.sampleName}/igv.all.splice {input.bam}
		"""


## Count k-mer frequency & Calculate k-mer correction scale factor
rule count_kmers:
	input:
		#sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz"
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/kmer.freq.txt"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		countKmersFromAlign.sh -l 3 -r 3 -g {genomeFa} -s {chrom_size} -f -v {input} > {output}
		"""

## Calculate k-mer adjustment scale factor based on genomewide frequency
rule calc_kmer_scale:
	input:
		sampleDir + "/{sampleName}/QC/kmer.freq.txt",
	output:
		sampleDir + "/{sampleName}/QC/kmer.scaleFactor.pseudo{kmer_pseudo}.txt"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		calcKmerScaleFactor.sh -g {kmer_genome} -p {kmer_pseudo} -v {input} > {output}
		"""


## 1bp-resolution bigwig file after k-mer correction
rule make_bigwig1bp_corrected:
	input:
		plus = sampleDir + "/{sampleName}/igv.1bp.plus.bw",
		minus = sampleDir + "/{sampleName}/igv.1bp.minus.bw",
		scale = sampleDir + "/{sampleName}/QC/kmer.scaleFactor.pseudo{kmer_pseudo}.txt",
		chrom = chrom_size
	output:
		sampleDir + "/{sampleName}/igv.1bp.corrected{kmer_pseudo}.plus.bw",
		sampleDir + "/{sampleName}/igv.1bp.corrected{kmer_pseudo}.minus.bw"
	message:
		"Making 1bp-resolution bigWig files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		correctKmerBiasBW.sh -l 3 -r 3 -g {genomeFa} -s {input.chrom} -k {input.scale} \
			-o {sampleDir}/{wildcards.sampleName}/igv.1bp.corrected{kmer_pseudo} -v \
			{sampleDir}/{wildcards.sampleName}/igv.1bp
		"""


## Create homer tag directories
rule make_tagdir_all:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		directory(sampleDir + "/{sampleName}/TSV.all"),
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.makeHomerDir.sh -o {output} -r 100 -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.frag}
		"""

rule make_tagdir_nfr:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		directory(sampleDir + "/{sampleName}/TSV.nfr"),
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.makeHomerDir.sh -o {output} -l 0 -L 119 -r 100 -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.frag}
		"""

rule make_tagdir_nuc:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		directory(sampleDir + "/{sampleName}/TSV.nuc")
	message:
		"Making Homer tag directory... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		ngs.makeHomerDir.sh -o {output} -l 151 -L 1000000 -r 100 -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.frag}
		"""


## Retrieve ctrl sample name for peak calling
def get_ctrl_name(sampleName):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	assert(len(ctrlName.tolist()[0]) != 0)
	return ctrlName.tolist()[0]


## Returns peak calling input tagDir(s): ctrl (optional) & target
def get_peakcall_input(sampleName, fragment, getCtrl=True):
	#ctrlName = samples.Ctrl[samples.Name == sampleName]
	#ctrlName = ctrlName.tolist()[0]
	if getCtrl:
		ctrlName = get_ctrl_name(sampleName)
	else:
		ctrlName = "NULL"

	if ctrlName.upper() == "NULL":
		return [ sampleDir + "/" + sampleName + "/TSV." + fragment ]
	else:
		return [ sampleDir + "/" + sampleName + "/TSV." + fragment, sampleDir + "/" + ctrlName + "/TSV." + fragment ]

## Return peak calling option from sample.tsv
def get_peakcall_opt(sampleName):
	#default options
	peakOpt = "-fragLength 100 -inputFragLength 100 -C 0"

	if "PeakOpt" in samples:
		optStr = samples.PeakOpt[samples.Name == sampleName]
		assert( len(optStr) == 1 )
		optStr = optStr.tolist()[0]
		if optStr != "NULL":
			peakOpt = peakOpt + " " + optStr

	return peakOpt


## Peak calling in factor mode using resized fragment
rule call_peaks_factor:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"nfr"),
		bw = sampleDir + "/{sampleName}/igv.nfr.con.bw",
	output:
		expand(sampleDir + "/{{sampleName}}/HomerPeak.factor/peak.exBL.1rpm.{ext}", ext=["bed", "stat"])
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.factor",
		optStr = lambda wildcards: "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {peak_mask} -b {input.bw} -s {params.optStr} {input.tagDir}
		"""

## Peak calling in factor mode using resized fragment
rule call_peaks_factor_no_ctrl:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"nfr", getCtrl=False),
		bw = sampleDir + "/{sampleName}/igv.nfr.con.bw",
	output:
		expand(sampleDir + "/{{sampleName}}/HomerPeak.factor.noCtrl/peak.exBL.1rpm.{ext}", ext=["bed", "stat"])
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.factor.noCtrl",
		optStr = lambda wildcards: "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {peak_mask} -b {input.bw} -s {params.optStr} {input.tagDir}
		"""


## Peak calling in factor mode using original size fragment
rule call_peaks_factor_allfrag:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"all"),
		bw = sampleDir + "/{sampleName}/igv.all.con.bw",
	output:
		expand(sampleDir + "/{{sampleName}}/HomerPeak.factor.allFrag/peak.exBL.1rpm.{ext}", ext=["bed","stat"])
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.factor.allFrag",
		optStr = lambda wildcards: "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {peak_mask} -b {input.bw} -s {params.optStr} {input.tagDir}
		"""

## Peak calling in histone mode using resized fragment
rule call_peaks_histone:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"nuc"),
	output:
		expand(sampleDir + "/{{sampleName}}/HomerPeak.histone/peak.exBL.{ext}", ext=["bed", "stat"])
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.histone",
		optStr = lambda wildcards: "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.peakCallHistone.sh -o {params.peakDir} -m {peak_mask} -s {params.optStr} {input.tagDir}
		"""

## Peak calling in histone mode using original size fragment
rule call_peaks_histone_allfrag:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"all")
	output:
		expand(sampleDir + "/{{sampleName}}/HomerPeak.histone.allFrag/peak.exBL.{ext}", ext=["bed","stat"])
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.histone.allFrag",
		optStr = lambda wildcards: "\"" + get_peakcall_opt(wildcards.sampleName) + "\""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.peakCallHistone.sh -o {params.peakDir} -m {peak_mask} -s {params.optStr} {input.tagDir}
		"""



rule run_homer_motif:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Motif/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/Motif/Homer.all
			touch {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/Motif/Homer.all/homerResults.html
			exit 0
		fi
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b {homerPreparseDir} \
			-o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/Motif/Homer.all {input}
		"""

rule run_homer_motif_allfrag:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Motif/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.allFrag/Motif/Homer.all
			touch {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.allFrag/Motif/Homer.all/homerResults.html
			exit 0
		fi
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b {homerPreparseDir} \
			-o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.allFrag/Motif/Homer.all {input}
		"""

rule run_meme_motif_rand5k:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
		#db = meme_db
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/Motif/MEME.random5k/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k peaks [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/Motif/MEME.random5k
			touch {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/Motif/MEME.random5k/meme-chip.html
			exit 0
		fi
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d {meme_db} \
			-o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/Motif/MEME.random5k {input.bed}
		"""


rule run_meme_motif_rand5k_allfrag:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed"
		#db = meme_db
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/Motif/MEME.random5k/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k peaks [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.allFrag/Motif/MEME.random5k
			touch {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.allFrag/Motif/MEME.random5k/meme-chip.html
			exit 0
		fi
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d {meme_db} \
			-o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.allFrag/Motif/MEME.random5k {input.bed}
		"""


#get bw file for control sample
#bwType variable should be either "nuc" or "nfr"
#this function assumes that a ctrl sample exists and is not NULL
def get_ctrl_bw(sampleName, bwType):
	ctrlName = get_ctrl_name(sampleName)
	return sampleDir + "/" + ctrlName + "/igv." + bwType + ".con.bw"

#get bw files for the sample and the corresponding ctrl
def get_bw_pairs(sampleName):
	
	#get bw files for sample and add them to list
	sampleNfrBW = sampleDir + "/" + sampleName + "/igv.nfr.con.bw"
	sampleNucBW = sampleDir + "/" + sampleName + "/igv.nuc.con.bw" 

	#get ctrl sample name
	ctrlName = get_ctrl_name(sampleName)
	
	#return sample bw only if ctrl = NULL
	if ctrlName == "NULL":
		return [sampleNfrBW, sampleNucBW]
	
	#return sample and ctrl bw 
	else:
		ctrlNfrBW = get_ctrl_bw(sampleName, "nfr")
		ctrlNucBW = get_ctrl_bw(sampleName, "nuc")
		return [sampleNfrBW, sampleNucBW, ctrlNfrBW, ctrlNucBW]


rule draw_peak_heatmap_factor:
	input:
		bed=sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/heatmap.exBL.1rpm.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0

		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		
		else

			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 2000 -b 20 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""


rule draw_peak_heatmap_histone:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/heatmap.exBL.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""

		module load Cutlery/1.0
		
		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		
		else

			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 10000 -b 100 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""


rule draw_peak_heatmap_factor_allfrag:
	input:
		bed=sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/heatmap.exBL.1rpm.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0

		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		
		else

			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 2000 -b 20 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""

rule draw_peak_heatmap_histone_allfrag:
	input:
		bed=sampleDir + "/{sampleName}/HomerPeak.histone.allFrag/peak.exBL.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone.allFrag/heatmap.exBL.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0

		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		else
			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 10000 -b 20 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""

rule analyze_footprint_homer:
	input:
		peak = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed",
		motif = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed.all.noBG/homerResults.html",
		bwPlus = sampleDir + "/{sampleName}/igv.1bp.plus.bw",
		bwMinus = sampleDir + "/{sampleName}/igv.1bp.minus.bw"
	output:
		sampleDir + "/{sampleName}/Footprint.Homer.default/cnr.4.complete"
	params:
		outPrefix = sampleDir + "/{sampleName}/Footprint.Homer.default/cnr",
		motifDir = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed.all.noBG",
		bwPrefix = sampleDir + "/{sampleName}/igv.1bp",
		cpu = cluster["analyze_footprint_homer"]["cpu"]
	message:
		"Analyzing CUT&RUN footprints... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.analyzeFootprintBatch.r -o {params.outPrefix} -n {wildcards.sampleName} -g {genome} -b {params.bwPrefix} -p {params.cpu} \
			{input.peak} {params.motifDir}
		"""


rule analyze_footprint_homer_corrected:
	input:
		peak = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed",
		motif = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed.all.noBG/homerResults.html",
		bwPlus = sampleDir + "/{sampleName}/igv.1bp.corrected{kmer_pseudo}.plus.bw",
		bwMinus = sampleDir + "/{sampleName}/igv.1bp.corrected{kmer_pseudo}.minus.bw"
	output:
		sampleDir + "/{sampleName}/Footprint.Homer.corrected{kmer_pseudo}/cnr.4.complete"
	params:
		outPrefix = sampleDir + "/{sampleName}/Footprint.Homer.corrected{kmer_pseudo}/cnr",
		motifDir = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed.all.noBG",
		bwPrefix = sampleDir + "/{sampleName}/igv.1bp.corrected{kmer_pseudo}",
		cpu = cluster["analyze_footprint_homer_corrected"]["cpu"]
	message:
		"Analyzing CUT&RUN footprints... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.analyzeFootprintBatch.r -o {params.outPrefix} -n {wildcards.sampleName} -g {genome} -b {params.bwPrefix} -p {params.cpu} \
			{input.peak} {params.motifDir}
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
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | cut -f 6`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		ngs.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.all} {input.all}
		ngs.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.nfr} {input.nfr}
		ngs.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.nuc} {input.nuc}
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
		module load Cutlery/1.0
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
		module load Cutlery/1.0
		scaleFactor=`cat {input.spikeinCnt} | gawk '$1=="'{wildcards.sampleName}'"' | gawk '{{ printf "%f", 100000/$3 }}'`
		if [ $scaleFactor == "" ];then
			echo -e "Error: empty scale factor" >&2
			exit 1
		fi
		ngs.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""

'''


#take sampleName and sample.tsv file as input
#return peakMode
def get_peak_mode(sampleName):
	peakMode = samples.PeakMode[samples.Name == sampleName]
	assert(len(peakMode.tolist()[0]) != 0)
	return peakMode.tolist()[0]

rule calc_frag_QC:
	input:
		fragDist = sampleDir + "/{sampleName}/QC/fragLen.dist.txt"
	output:
		sampleDir + "/{sampleName}/QC/fragMix.txt"
	message:
		"Performing fragment QC on sample... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.calcFragMixture.py -o {sampleDir}/{wildcards.sampleName}/QC/fragMix \
		{input.fragDist}
		"""

#take sampleName as input
#return group
def get_group(sampleName):
	group = samples.Group[samples.Name == sampleName]
	assert(len(group.tolist()[0]) != 0)
	return group.tolist()[0]

## measure % of fragments mapped to promoter regions
rule measure_promoter_fraction:
	input:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/promoter_portion.txt"
	message:
		"Checking fragments within promoters... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		N_all=`zcat {input} | wc -l`
		N_prom=`zcat {input} \\
			| gawk 'BEGIN{{ n=0 }}{{
					c=($2+$3)/2
					if(c<50) c=50
					printf "%s\\t%d\\t%d\\n", $1,c-50,c+50
				}}' \\
			| intersectBed -a stdin -b {bed_promoter} -u \\
			| wc -l`
		echo -e "${{N_all}}\\t${{N_prom}}" | gawk '{{ printf "Name\\t{wildcards.sampleName}\\nTTC\\t%d\\nPromoter\\t%d\\nPercentage\\t%.1f\\n", $1,$2,$2/$1*100 }}' > {output}
		"""


rule make_promotercnt_table:
	input:
		expand(sampleDir + "/{sampleName}/QC/promoter_portion.txt", sampleName=samples.Name.tolist())
	output:
		expand(qcDir + "/promoterCnt.{ext}", ext=[ "txt", "pdf", "png" ])
	message:
		"Making promoter count table..."
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.makePromCntTable.r -o {qcDir}/promoterCnt {input}
		"""



################ MACS Peak Calling in PE mode #######################
## Return bam file (with full path) for a given sample name
# Return control bam file if mode == "ctrl"
# NOTE: needs to handle no-ctrl case

rule make_bam_nfr:
	input:
		lambda wildcards: get_bam_for_downstream(wildcards.sampleName)
	output:
		( dedupDir if doDedup else alignDir ) + "/{sampleName}/Split/align.nfr.bam"
	message:
		"Making NFR bam file... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load samtools/1.14.0
		samtools view -h {input} \
			| gawk 'substr($0,1,1)=="@" || ($9 < 120 && $9 > -120)' \
			| samtools view -b \
			> {output}
		"""

rule make_bam_nuc:
	input:
		lambda wildcards: get_bam_for_downstream(wildcards.sampleName)
	output:
		( dedupDir if doDedup else alignDir ) + "/{sampleName}/Split/align.nuc.bam"

	message:
		"Making NUC bam file... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load samtools/1.14.0
		samtools view -h {input} \
			| gawk '/^@/ || ($9 >= 120 || $9 <= -120)' \
			| samtools view -b \
			> {output}
		"""

## find control sample name for peak calling using target sample name
def get_ctrl_name(sampleName):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	ctrlName = ctrlName.tolist()[0]
	return ctrlName


def get_bam_for_macs(sampleName, fragment, mode="target"):
	assert fragment in [ "nfr", "nuc", "all" ]
	assert mode in [ "target", "ctrl" ]
	bamDir = dedupDir if doDedup else alignDir
	if mode == "target":
		name = sampleName
	else:
		name = get_ctrl_name(sampleName)
	
	if fragment in [ "nfr", "nuc" ]:
		bam = bamDir + "/" + name + "/Split/align." + fragment + ".bam"
	else:
		bam = bamDir + "/" + name + "/align.bam"
	return bam

def get_frag_for_macs(sampleName, fragment, mode="target"):
	assert fragment in [ "nfr", "nuc", "all" ]
	assert mode in [ "target", "ctrl" ]
	if mode == "target":
		name = sampleName
	else:
		name = get_ctrl_name(sampleName)
	
	if fragment in [ "nfr", "nuc" ]:
		frag = sampleDir + "/" + name + "/fragment." + fragment + ".bed.gz"
	else:
		frag = sampleDir + "/" + name + "/fragment.bed.gz"
	return frag

## MACS peak calling vs control: NFR
rule call_peak_macs_factor:
	input:
		target = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "nfr", mode="target"),
		ctrl = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "nfr", mode="ctrl")
	output:
		peak = sampleDir + "/{sampleName}/MACS2.factor/{sampleName}_summits.exBL.bed",
		log = sampleDir + "/{sampleName}/MACS2.factor/{sampleName}.log"
	message:
		"Calling TF peaks using MACS.. [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MACS/2.2.9.1
		module load bedtools/2.27.0
		macs2 callpeak -t {input.target} -c {input.ctrl} -f BAMPE -n {wildcards.sampleName} --outdir {params.outDir} -g {species_macs} --keep-dup all --call-summits 2>&1 | tee {output.log}
		intersectBed -a {params.outDir}/{wildcards.sampleName}_summits.bed -b {params.mask} -v \
			| gawk '{{ printf "%s\\t%d\\t%d\\t%s\\t%s\\t+\\n", $1,$2,$3,$4,$5 }}' \
			| sort -k5,5nr \
			> {output.peak}
		"""

rule draw_peak_heatmap_factor_macs:
	input:
		bed = sampleDir + "/{sampleName}/MACS2.factor/{sampleName}_summits.exBL.bed",
		bw = lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/MACS2.factor/heatmap.exBL.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0

		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		
		else

			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 2000 -b 20 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""

# MACS peak calling vs control: all fragment
rule call_peak_macs_factor_allfrag:
	input:
		target = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "all", mode="target"),
		ctrl = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "all", mode="ctrl")
	output:
		peak = sampleDir + "/{sampleName}/MACS2.factor.allFrag/{sampleName}_summits.exBL.bed",
		log = sampleDir + "/{sampleName}/MACS2.factor.allFrag/{sampleName}.log"
	message:
		"Calling TF peaks using MACS.. [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MACS/2.2.9.1
		module load bedtools/2.27.0
		macs2 callpeak -t {input.target} -c {input.ctrl} -f BAMPE -n {wildcards.sampleName} --outdir {params.outDir} -g {species_macs} --keep-dup all --call-summits 2>&1 | tee {output.log}
		intersectBed -a {params.outDir}/{wildcards.sampleName}_summits.bed -b {params.mask} -v \
			| gawk '{{ printf "%s\\t%d\\t%d\\t%s\\t%s\\t+\\n", $1,$2,$3,$4,$5 }}' \
			| sort -k5,5nr \
			> {output.peak}
		"""

rule draw_peak_heatmap_factor_macs_allfrag:
	input:
		bed = sampleDir + "/{sampleName}/MACS2.factor.allFrag/{sampleName}_summits.exBL.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/MACS2.factor.allFrag/heatmap.exBL.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0

		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		
		else

			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 2000 -b 20 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""


## MACS peak calling vs control: NUC
rule call_peak_macs_histone:
	input:
		target = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "nuc", mode="target"),
		ctrl = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "nuc", mode="ctrl")
	output:
		peak = sampleDir + "/{sampleName}/MACS2.histone/{sampleName}_broad.exBL.bed",
		log = sampleDir + "/{sampleName}/MACS2.histone/{sampleName}.log"
	message:
		"Calling histone peaks using MACS.. [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MACS/2.2.9.1
		module load bedtools/2.27.0
		macs2 callpeak -t {input.target} -c {input.ctrl} -f BAMPE -n {wildcards.sampleName} --outdir {params.outDir} -g {species_macs} --broad --keep-dup all 2>&1 | tee {output.log}
		
		grep -v '^#' {sampleDir}/{wildcards.sampleName}/MACS2.histone/{wildcards.sampleName}_peaks.xls \
			| tail -n +3 \
			| intersectBed -a stdin -b {params.mask} -v \
			| gawk '{{ printf "%s\\t%d\\t%d\\t%s\\t%d\\t+\\n", $1,$2,$3,$9,$4 }}' \
			| sort -k5,5nr \
			> {output.peak}
		"""

rule draw_peak_heatmap_histone_macs:
	input:
		bed = sampleDir + "/{sampleName}/MACS2.histone/{sampleName}_broad.exBL.bed",
		bw = lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/MACS2.histone/heatmap.exBL.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0

		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		
		else

			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 10000 -b 20 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""

## MACS histone peak calling vs control: all fragment
rule call_peak_macs_histone_allfrag:
	input:
		target = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "all", mode="target"),
		ctrl = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "all", mode="ctrl")
	output:
		peak = sampleDir + "/{sampleName}/MACS2.histone.allFrag/{sampleName}_broad.exBL.bed",
		log = sampleDir + "/{sampleName}/MACS2.histone.allFrag/{sampleName}.log"
	message:
		"Calling histone peaks using MACS.. [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MACS/2.2.9.1
		module load bedtools/2.27.0
		macs2 callpeak -t {input.target} -c {input.ctrl} -f BAMPE -n {wildcards.sampleName} --outdir {params.outDir} -g {species_macs} --broad --keep-dup all 2>&1 | tee {output.log}

		grep -v '^#' {sampleDir}/{wildcards.sampleName}/MACS2.histone.allFrag/{wildcards.sampleName}_peaks.xls \
			| tail -n +3 \
			| intersectBed -a stdin -b {params.mask} -v \
			| gawk '{{ printf "%s\\t%d\\t%d\\t%s\\t%d\\t+\\n", $1,$2,$3,$9,$4 }}' \
			| sort -k5,5nr \
			> {output.peak}
		"""

rule draw_peak_heatmap_histone_macs_allfrag:
	input:
		bed = sampleDir + "/{sampleName}/MACS2.histone.allFrag/{sampleName}_broad.exBL.bed",
		bw = lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/MACS2.histone.allFrag/heatmap.exBL.png"
	params:
		outPrefix = lambda wildcards, output: __import__("re").sub(".png$","", output[0])
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0

		## Check if the file is empty; if empty, create image saying No peak detected
		if [ ! -s "{input.bed}" ]; then
			convert -size 800x600 xc:white -gravity Center -pointsize 48 -font ${{CUTLERY}}/Resource/freeserif/FreeSerif.otf -annotate 0 "No peak detected" {output}
		
		else

			cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 10000 -b 20 \
				-o {params.outPrefix} \
				{input.bed} {input.bw}
		fi
		"""

## MACS peak calling with relaxed criteria using p-value for IDR analysis
rule call_peak_macs_factor_relax:
	input:
		target = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "nfr", mode="target"),
		ctrl = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "nfr", mode="ctrl")
	output:
		peak = sampleDir + "/{sampleName}/MACS2.factor.relax/{sampleName}_summits.exBL.bed",
		peak2 = sampleDir + "/{sampleName}/MACS2.factor.relax/{sampleName}_peaks.sorted.narrowPeak",
		log = sampleDir + "/{sampleName}/MACS2.factor.relax/{sampleName}.log"
	message:
		"Calling TF peaks using MACS.. [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MACS/2.2.9.1
		module load bedtools/2.27.0
		macs2 callpeak -t {input.target} -c {input.ctrl} -f BAMPE -n {wildcards.sampleName} \
			--outdir {params.outDir} -g {species_macs} --keep-dup all --call-summits -p 0.001 \
			2>&1 | tee {output.log}
		intersectBed -a {params.outDir}/{wildcards.sampleName}_summits.bed -b {params.mask} -v > {output.peak}
		sort -k8,8nr {params.outDir}/{wildcards.sampleName}_peaks.narrowPeak > {output.peak2}
		"""

## MACS peak calling without control
rule call_peak_macs_factor_wo_ctrl:
	input:
		target = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "nfr", mode="target")
	output:
		peak = sampleDir + "/{sampleName}/MACS2.factor.wo_ctrl/{sampleName}_summits.exBL.bed",
		log = sampleDir + "/{sampleName}/MACS2.factor.wo_ctrl/{sampleName}.log"
	message:
		"Calling TF peaks/SE ... [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MACS/2.2.9.1
		module load bedtools/2.27.0
		macs2 callpeak -t {input.target} -f BAMPE -n {wildcards.sampleName} --outdir {params.outDir} -g {species_macs} --keep-dup all --call-summits 2>&1 | tee {output.log}
		intersectBed -a {params.outDir}/{wildcards.sampleName}_summits.bed -b {params.mask} -v > {output.peak}
		"""

rule call_peak_macs_factor_allfrag_wo_ctrl:
	input:
		target = lambda wildcards: get_bam_for_macs(wildcards.sampleName, fragment = "all", mode="target")
	output:
		peak = sampleDir + "/{sampleName}/MACS2.factor.allFrag.wo_ctrl/{sampleName}_summits.exBL.bed",
		log = sampleDir + "/{sampleName}/MACS2.factor.allFrag.wo_ctrl/{sampleName}.log"
	message:
		"Calling TF peaks/SE ... [{wildcards.sampleName}]"
	params:
		mask = peak_mask,
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MACS/2.2.9.1
		module load bedtools/2.27.0
		macs2 callpeak -t {input.target} -f BAMPE -n {wildcards.sampleName} --outdir {params.outDir} -g {species_macs} --keep-dup all --call-summits 2>&1 | tee {output.log}
		intersectBed -a {params.outDir}/{wildcards.sampleName}_summits.bed -b {params.mask} -v > {output.peak}
		"""


rule macs_run_homer_motif:
	input:
		sampleDir + "/{sampleName}/MACS2.factor/{sampleName}_summits.exBL.bed",
	output:
		sampleDir + "/{sampleName}/MACS2.factor/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load Motif/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {params.outDir}
			touch {params.outDir}/homerResults.html
			exit 0
		fi
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b {homerPreparseDir} \
			-o {params.outDir} {input}
		"""

rule macs_run_homer_motif_allFrag:
	input:
		sampleDir + "/{sampleName}/MACS2.factor.allFrag/{sampleName}_summits.exBL.bed",
	output:
		sampleDir + "/{sampleName}/MACS2.factor.allFrag/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load Motif/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {params.outDir}
			touch {params.outDir}/homerResults.html
			exit 0
		fi
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b {homerPreparseDir} \
			-o {params.outDir} {input}
		"""


rule macs_run_meme_motif_rand5k:
	input:
		sampleDir + "/{sampleName}/MACS2.factor/{sampleName}_summits.exBL.bed"
	output:
		sampleDir + "/{sampleName}/MACS2.factor/Motif/MEME.random5k/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k peaks [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {params.outDir}
			touch {params.outDir}/meme-chip.html
			exit 0
		fi
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d {meme_db} \
			-o {params.outDir} {input}
		"""

rule macs_run_meme_motif_rand5k_allFrag:
	input:
		sampleDir + "/{sampleName}/MACS2.factor.allFrag/{sampleName}_summits.exBL.bed"
	output:
		sampleDir + "/{sampleName}/MACS2.factor.allFrag/Motif/MEME.random5k/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k peaks [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {params.outDir}
			touch {params.outDir}/meme-chip.html
			exit 0
		fi
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d {meme_db} \
			-o {params.outDir} {input}
		"""



rule macs_run_homer_motif_wo_ctrl:
	input:
		sampleDir + "/{sampleName}/MACS2.factor.wo_ctrl/{sampleName}_summits.exBL.bed",
	output:
		sampleDir + "/{sampleName}/MACS2.factor.wo_ctrl/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load Motif/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {params.outDir}
			touch {params.outDir}/homerResults.html
			exit 0
		fi
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b {homerPreparseDir} \
			-o {params.outDir} {input}
		"""


rule macs_run_homer_motif_allfrag_wo_ctrl:
	input:
		sampleDir + "/{sampleName}/MACS2.factor.allFrag.wo_ctrl/{sampleName}_summits.exBL.bed",
	output:
		sampleDir + "/{sampleName}/MACS2.factor.allFrag.wo_ctrl/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	params:
		outDir = lambda wildcards, output: __import__("os").path.dirname(output[0])
	shell:
		"""
		module purge
		module load Motif/1.0
		n=`cat {input} | wc -l`
		if [ $n -eq 0 ];then
			mkdir -p {params.outDir}
			touch {params.outDir}/homerResults.html
			exit 0
		fi
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b {homerPreparseDir} \
			-o {params.outDir} {input}
		"""


########################################
## Rules for report generation
########################################

#take sampleName as input
#return path to peak examples and heatmap png
def get_peak_example(sampleName):
	peakMode = get_peak_mode(sampleName)
	if peakMode == "histone":
		return sampleDir + "/" + sampleName + "/HomerPeak.histone/peak.examples.png"
	elif peakMode == "factor":
		return sampleDir + "/" + sampleName + "/HomerPeak.factor/peak.examples.png"
	elif peakMode == "NULL":
		return ""

def get_heatmap(sampleName):
	peakMode = get_peak_mode(sampleName)
	if peakMode == "histone":
		return sampleDir + "/" + sampleName + "/HomerPeak.histone/heatmap.exBL.png"
	elif peakMode == "factor":
		return sampleDir + "/" + sampleName + "/HomerPeak.factor/heatmap.exBL.1rpm.png"
	elif peakMode == "NULL":
		return ""

def get_peakIntersectPerc(sampleName):
	peakMode = get_peak_mode(sampleName)
	if peakMode == "histone":
		return sampleDir + "/" + sampleName + "/HomerPeak.histone/peak.exBL.stat"
	elif peakMode == "factor":
		return sampleDir + "/" + sampleName + "/HomerPeak.factor/peak.exBL.1rpm.stat"
	elif peakMode == "NULL":
		return ""

rule create_report_per_sample:
	input:
		alnStat = qcDir + "/alignStat.txt",
		uniqFragCnt = qcDir + "/uniqFragCnt.txt",
		baseFreqPNG = sampleDir + "/{sampleName}/QC/base_freq.png",
		fragDistPNG = sampleDir + "/{sampleName}/QC/fragLen.dist.png",
		fragQC = sampleDir + "/{sampleName}/QC/fragMix.txt",
		heatmap = lambda wildcards: get_heatmap(wildcards.sampleName),
		peakStat = lambda wildcards: get_peakIntersectPerc(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/QC/Report.html"
	message:
		"Creating report for sample... [{wildcards.sampleName}]"
	params:
		group = lambda wildcards: get_group(wildcards.sampleName)
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.createSampleReportHTML.r -o Report -t {src_sampleInfo} -g {params.group} -s {sampleDir}/{wildcards.sampleName} -q {qcDir}
		"""


rule create_final_report:
	input:
		uniqFragCnt = qcDir + "/uniqFragCnt.txt",
		fragDist = expand(sampleDir + "/{sampleName}/QC/fragLen.dist.txt", sampleName=samples.Name.tolist()),
		fragQC = expand(sampleDir + "/{sampleName}/QC/fragMix.txt", sampleName=samples.Name.tolist()),
		histoneHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.histone/heatmap.exBL.png", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		factorHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.factor/heatmap.exBL.1rpm.png", sampleName = samples.Name[samples.PeakMode=="factor"].tolist()),
		histonePeakIntersectPerc = expand(sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.stat", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		factorPeakIntersectPerc = expand(sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.stat", sampleName = samples.Name[samples.PeakMode=="factor"].tolist())
	output:
		"Report.html"
	message:
		"Creating final report in HTML..."
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.createReportHTML.r -o Report -t {src_sampleInfo} -s {sampleDir} -q {qcDir}
		"""

rule create_report_per_sample_pooled:
	input:
		uniqFragCnt = qcDir + "/uniqFragCnt.txt",
		baseFreqPNG = sampleDir + "/{sampleName}/QC/base_freq.png",
		fragDistPNG = sampleDir + "/{sampleName}/QC/fragLen.dist.png",
		fragQC = sampleDir + "/{sampleName}/QC/fragMix.txt",
		heatmap = lambda wildcards: get_heatmap(wildcards.sampleName),
		peakStat = lambda wildcards: get_peakIntersectPerc(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/QC/Report_pooled.html"
	message:
		"Creating report for sample... [{wildcards.sampleName}]"
	params:
		group = lambda wildcards: get_group(wildcards.sampleName)
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.createPooledSampleReportHTML.r -o Report_pooled -t {src_sampleInfo} -g {params.group} -s {sampleDir}/{wildcards.sampleName} -q {qcDir}		
		"""

rule create_final_report_pooled:
	input:
		fragDist = expand(sampleDir + "/{sampleName}/QC/fragLen.dist.txt", sampleName=samples.Name.tolist()),
		fragQC = expand(sampleDir + "/{sampleName}/QC/fragMix.txt", sampleName=samples.Name.tolist()),
		histoneHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.histone/heatmap.exBL.png", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		factorHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.factor/heatmap.exBL.1rpm.png", sampleName = samples.Name[samples.PeakMode=="factor"].tolist()),
		histonePeakIntersectPerc = expand(sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.stat", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		factorPeakIntersectPerc = expand(sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.stat", sampleName = samples.Name[samples.PeakMode=="factor"].tolist())
	output:
		"Report_pooled.html"
	message:
		"Creating final report in HTML..."
	shell:
		"""
		module purge
		module load Cutlery/1.0
		cnr.createPooledReportHTML.r -o Report_pooled -t {src_sampleInfo} -s {sampleDir} -q {qcDir} -f fragMix
		"""
