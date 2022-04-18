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


if "meme_db" not in locals():
	meme_db = os.environ["LIMLAB_BASE"] + "/Motif/MEME_DB/Merged_By_Lim.meme"

if "numHighestPeaks" not in locals():
	numHighestPeaks = 5


#################################
## Rule Start

## bamDir is used to designated bam file directory for downstream analysis
## Originally motivated for reuse or rule.post.smk in replicate pooling in Snakemake.pool
## If bamDir is not defined, i.e. for simply processing individual replicates not pooling
## dedupDir or alignDir is selected 
if "bamDir" not in locals():
	if doDedup:
		bamDir = dedupDir
	else:
		bamDir = alignDir

## Convert BAM to fragment bed file
## - Chromosome filtering
## - 0x2  : Concordant pairs only
## - 0x400: Removes duplicates
rule make_fragment:
	input:
		bam = bamDir + "/{sampleName}/align.bam",
		bai = bamDir + "/{sampleName}/align.bam.bai"
	output:
		sampleDir + "/{sampleName}/fragment.bed.gz"
	message:
		"Making fragment bed files... [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load Cutlery/1.0
		mkdir -p {sampleDir}/{wildcards.sampleName}
		ngs.bamToFragment.py -c "{chrRegexAll}" -f 0x2 -F 0x400 {input.bam} | sort -S 2G -k1,1 -k2,2n -k3,3n | gzip > {output}
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
		module load Cutlery/1.0
		ngs.makeUniqCntTable.r -o {qcDir}/uniqFragCnt {input}
		"""

## Nucleotide base frequence around 5'-ends of read 1 & 2
## NEED REVISION for OUTPUT names & directories
rule check_baseFreq:
	input:
		#frag = sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz",
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		genomeFa = genomeFa
		#bamDir + "/{sampleName}.bam"
	output:
		sampleDir + "/{sampleName}/QC/base_freq.png",
		sampleDir + "/{sampleName}/QC/base_freq.html"
#		sampleDir + "/{sampleName}/QC/base_freq.R1.png",
#		sampleDir + "/{sampleName}/QC/base_freq.R2.png"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		checkBaseFreq.plot.sh -o {sampleDir}/{wildcards.sampleName}/QC/base_freq \
			-n {wildcards.sampleName} -g {input.genomeFa} -c "{chrRegexTarget}" -m both -l 20 -f -i -v {input.frag}
		"""
#		bamToBed.separate.sh -o {sampleDir}/{wildcards.sampleName}/tmp {input}
#		checkBaseFreq.plot.sh -g {genomeFa} -n {wildcards.sampleName} -o {sampleDir}/{wildcards.sampleName}/QC/base_freq.R1 {sampleDir}/{wildcards.sampleName}/tmp.R1.bed.gz
#		checkBaseFreq.plot.sh -g {genomeFa} -n {wildcards.sampleName} -o {sampleDir}/{wildcards.sampleName}/QC/base_freq.R2 {sampleDir}/{wildcards.sampleName}/tmp.R2.bed.gz
#		rm {sampleDir}/{wildcards.sampleName}/tmp.R1.bed.gz
#		rm {sampleDir}/{wildcards.sampleName}/tmp.R2.bed.gz


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
		module load Cutlery/1.0
		cnr.splitBamToBed.sh -o {sampleDir}/{wildcards.sampleName}/Fragments/frag -c "{chrRegexTarget}" {input}
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
		module load Cutlery/1.0
		cnr.drawAutoCorFrag.r -o {sampleDir}/{wildcards.sampleName}/QC/enrich {input}
		"""

rule count_spikein:
	input:
		#bamDir + "/{sampleName}.bam"
		sampleDir + "/{sampleName}/fragment.bed.gz"
	output:
		sampleDir + "/{sampleName}/QC/spikeCnt.txt"
	message:
		"Counting spikein tags... [{wildcards.sampleName}]"
	shell:
		"""
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
		module load Cutlery/1.0
		ngs.makeSpikeCntTable.r -o {qcDir}/spikein {input}
		"""

## NOTE: input fragments are already chromosome-filtered in split step
##	However, planning to add chrRegex option to cnr.fragToBigWig for better readibility & compatibility
rule make_bigwig:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		#all = sampleDir + "/{sampleName}/Fragments/frag.all.ctr.bed.gz"
		#nfr = sampleDir + "/{sampleName}/Fragments/frag.nfr.ctr.bed.gz",
		#nuc = sampleDir + "/{sampleName}/Fragments/frag.nuc.ctr.bed.gz",
		chrom = chrom_size
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
		module load Cutlery/1.0
		cnr.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -r 100 -m 5G -o {output.all} {input.frag}
		cnr.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 0 -L 119 -r 100 -m 5G -o {output.nfr} {input.frag}
		cnr.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 151 -L 1000000 -r 100 -m 5G -o {output.nuc} {input.frag}
		"""


## 1bp-resolution bigwig files
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
		module load Cutlery/1.0
		cnr.fragToBigWigStranded1bp.sh -o {sampleDir}/{wildcards.sampleName}/igv.1bp -g {input.chrom} -c "{chrRegexTarget}" -s 0 -m 5G {input.frag}
		"""
#		ngs.alignToBigWig.sh -o {sampleDir}/{wildcards.sampleName}/igv.1bp -g {chrom_size} -l 1 -m 5G -c "{chrRegexTarget}" {input}


## Raw read-scale / nonnegative tracks (positive values even for minus track)
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
		module load Cutlery/1.0
		cnr.fragToBigWigStranded1bp.sh -o {sampleDir}/{wildcards.sampleName}/igv.1bp.raw.abs -g {input.chrom} -c "{chrRegexTarget}" -s 1 -m 5G -n {input.frag}
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
		module load Cutlery/1.0
		countKmersFromAlign.sh -l 3 -r 3 -g {genomeFa} -s {chrom_size} -f -v {input} > {output}
		"""

rule calc_kmer_scale:
	input:
		sampleDir + "/{sampleName}/QC/kmer.freq.txt",
	output:
		sampleDir + "/{sampleName}/QC/kmer.scaleFactor.pseudo{kmer_pseudo}.txt"
	message:
		"Checking baseFrequency... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		calcKmerScaleFactor.sh -g {kmer_genome} -p {kmer_pseudo} -v {input} > {output}
		"""


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
		module load Cutlery/1.0
		correctKmerBiasBW.sh -l 3 -r 3 -g {genomeFa} -s {input.chrom} -k {input.scale} \
			-o {sampleDir}/{wildcards.sampleName}/igv.1bp.corrected{kmer_pseudo} -v \
			{sampleDir}/{wildcards.sampleName}/igv.1bp
		"""


rule make_bigwig_allfrag:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz",
		#all=sampleDir + "/{sampleName}/Fragments/frag.all.con.bed.gz"
		#nfr=sampleDir + "/{sampleName}/Fragments/frag.nfr.con.bed.gz",
		#nuc=sampleDir + "/{sampleName}/Fragments/frag.nuc.con.bed.gz",
		chrom = chrom_size
	output:
		all=sampleDir + "/{sampleName}/igv.all.con.bw",
		nfr=sampleDir + "/{sampleName}/igv.nfr.con.bw",
		nuc=sampleDir + "/{sampleName}/igv.nuc.con.bw"
	message:
		"Making bigWig files... [{wildcards.sampleName}]"
#	params:
#		memory = "5G"
	shell:
		"""
		module load Cutlery/1.0
		cnr.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -m 5G -o {output.all} {input.frag}
		cnr.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 0 -L 119 -m 5G -o {output.nfr} {input.frag}
		cnr.fragToBigWig.sh -g {input.chrom} -c "{chrRegexTarget}" -l 151 -L 1000000 -m 5G -o {output.nuc} {input.frag}
		"""

rule make_tagdir:
	input:
		frag = sampleDir + "/{sampleName}/fragment.bed.gz"
		#all = sampleDir + "/{sampleName}/Fragments/frag.all.ctr.bed.gz"
		#nfr = sampleDir + "/{sampleName}/Fragments/frag.nfr.ctr.bed.gz",
		#nuc = sampleDir + "/{sampleName}/Fragments/frag.nuc.ctr.bed.gz"
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
		module load Cutlery/1.0
		cnr.makeHomerDir.sh -o {output.all} -r 100 -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.frag}
		cnr.makeHomerDir.sh -o {output.nfr} -l 0 -L 119 -r 100 -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.frag}
		cnr.makeHomerDir.sh -o {output.nuc} -l 151 -L 1000000 -r 100 -n {wildcards.sampleName} -c "{chrRegexTarget}" {input.frag}
		"""


def get_ctrl_name(sampleName):
	ctrlName = samples.Ctrl[samples.Name == sampleName]
	assert(len(ctrlName.tolist()[0]) != 0)
	return ctrlName.tolist()[0]


## Returns peak calling input tagDir(s): ctrl (optional) & target
def get_peakcall_input(sampleName, fragment):
	#ctrlName = samples.Ctrl[samples.Name == sampleName]
	#ctrlName = ctrlName.tolist()[0]
	ctrlName = get_ctrl_name(sampleName)
	if ctrlName.upper() == "NULL":
		return [ sampleDir + "/" + sampleName + "/TSV." + fragment ]
	else:
		return [ sampleDir + "/" + ctrlName + "/TSV." + fragment, sampleDir + "/" + sampleName + "/TSV." + fragment ]


rule call_peaks_factor:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"nfr"),
		bw = sampleDir + "/{sampleName}/igv.nfr.con.bw",
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.factor",
		optStr = lambda wildcards, input: "-i" if len(input.tagDir)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {peak_mask} -b {input.bw} -s \"-fragLength 100 -inputFragLength 100\" {params.optStr} {input.tagDir}
		"""


rule call_peaks_factor_allfrag:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"all"),
		bw = sampleDir + "/{sampleName}/igv.all.con.bw",
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed"
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.factor.allFrag",
		optStr = lambda wildcards, input: "-i" if len(input.tagDir)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.peakCallTF.sh -o {params.peakDir} -m {peak_mask} -b {input.bw} -s \"-fragLength 100 -inputFragLength 100\" {params.optStr} {input.tagDir}
		"""


rule call_peaks_histone:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"nuc"),
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed"
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.histone",
		optStr = lambda wildcards, input: "-i" if len(input.tagDir)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.peakCallHistone.sh -o {params.peakDir} -m {peak_mask} -s \"-fragLength 100 -inputFragLength 100 -C 0\" {params.optStr} {input.tagDir}
		"""


rule call_peaks_histone_allfrag:
	input:
		tagDir = lambda wildcards: get_peakcall_input(wildcards.sampleName,"all")
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone.allFrag/peak.exBL.bed"
	params:
		peakDir = sampleDir + "/{sampleName}/HomerPeak.histone.allFrag",
		optStr = lambda wildcards, input: "-i" if len(input.tagDir)>1 else ""
	message:
		"Peak calling using Homer... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.peakCallHistone.sh -o {params.peakDir} -m {peak_mask} -s \"-fragLength 100 -inputFragLength 100 -C 0\" {params.optStr} {input.tagDir}
		"""


'''
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
rule run_homermotif_allfrag:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed.all.noBG/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load Motif/1.0
		runHomerMotif.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse -o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor.allFrag {input}
		"""
'''


rule run_homer_motif:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/Motif/Homer.all/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load Motif/1.0
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse \
			-o {sampleDir}/{wildcards.sampleName}/Motif/Homer.all {input}
		"""

rule run_homermotif_allfrag:
	input:
		sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed"
	output:
		sampleDir + "/{sampleName}/Motif/Homer.all.allFrag/homerResults.html"
	message:
		"Running Homer motif search... [{wildcards.sampleName}]"
	shell:
		"""
		module load Motif/1.0
		runHomerMotifSingle.sh -g {genome} -s 200 -p 4 -b /data/limlab/Resource/Homer.preparse \
			-o {sampleDir}/{wildcards.sampleName}/Motif/Homer.all.allFrag {input}
		"""

rule run_meme_motif_rand5k:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed"
		#db = meme_db
	output:
		sampleDir + "/{sampleName}/Motif/MEME.random5k/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k TSS peaks [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d {meme_db} \
			-o {sampleDir}/{wildcards.sampleName}/Motif/MEME.random5k {input.bed}
		"""


rule run_meme_motif_rand5k_allfrag:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.factor.allFrag/peak.exBL.1rpm.bed"
		#db = meme_db
	output:
		sampleDir + "/{sampleName}/Motif/MEME.random5k.allFrag/meme-chip.html"
	message:
		"Running MEME-ChIP motif search for random 5k TSS peaks [{wildcards.sampleName}]"
	shell:
		"""
		module purge
		module load MotifMEME/1.0
		runMemeChipSingle.sh -g {genomeFa} -s 200 -p 4 -r 5000 -d {meme_db}} \
			-o {sampleDir}/{wildcards.sampleName}/Motif/MEME.random5k.allFrag {input.bed}
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
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 2000 -b 20 \
			-o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/heatmap.exBL.1rpm \
			{input.bed} {input.bw}
		"""


rule draw_peak_heatmap_histone:
	input:
		bed = sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/HomerPeak.histone/heatmap.exBL.png"
	message:
		"Drawing peak profile heatmap... [{wildcards.sampleName}]"
	shell:
		"""
		module load Cutlery/1.0
		cnr.drawPeakHeatmap.r -t {wildcards.sampleName} -w 10000 -b 100 \
			-o {sampleDir}/{wildcards.sampleName}/HomerPeak.histone/heatmap.exBL \
			{input.bed} {input.bw}
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
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.all} {input.all}
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.nfr} {input.nfr}
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output.nuc} {input.nuc}
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
		cnr.fragToBigWig.sh -g {chrom_size} -m 5G -s $scaleFactor -o {output} {input.bed}
		"""

'''


#take sampleName and sample.tsv file as input
#return peakMode
def get_peak_mode(sampleName):
	peakMode = samples.PeakMode[samples.Name == sampleName]
	assert(len(peakMode.tolist()[0]) != 0)
	return peakMode.tolist()[0]


rule draw_peak_examples_histone:
	input:
		peak = sampleDir + "/{sampleName}/HomerPeak.histone/peak.exBL.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		expand(sampleDir + "/{{sampleName}}/HomerPeak.histone/peak.examples.{ext}", ext=["png", "pdf"])
	message:
		"Getting snapshot of highest scoring peaks... [{wildcards.sampleName}]"
	params:
		peakMode = lambda wildcards: get_peak_mode(wildcards.sampleName)
	shell:
		"""
		module load Cutlery/1.0
		cnr.drawPeakExamples.r -o {sampleDir}/{wildcards.sampleName}/HomerPeak.histone/peak.examples \
			-m {params.peakMode} -n "{numHighestPeaks}" -p {input.peak} {input.bw}
		"""


rule draw_peak_examples_factor:
	input:
		peak = sampleDir + "/{sampleName}/HomerPeak.factor/peak.exBL.1rpm.bed",
		bw=lambda wildcards: get_bw_pairs(wildcards.sampleName)
	output:
		expand(sampleDir + "/{{sampleName}}/HomerPeak.factor/peak.examples.{ext}", ext=["png", "pdf"])
	message:
		"Getting snapshot of highest scoring peaks... [{wildcards.sampleName}]"
	params:
		peakMode = lambda wildcards: get_peak_mode(wildcards.sampleName)
	shell:
		"""
		module load Cutlery/1.0
		cnr.drawPeakExamples.r -o {sampleDir}/{wildcards.sampleName}/HomerPeak.factor/peak.examples \
			-m {params.peakMode} -n "{numHighestPeaks}" -p {input.peak} {input.bw}
		"""

rule calc_frag_QC:
	input:
		fragDist = sampleDir + "/{sampleName}/QC/fragLen.dist.txt"
	output:
		sampleDir + "/{sampleName}/QC/fragMix.txt"
	message:
		"Performing fragment QC on sample... [{wildcards.sampleName}]"
	shell:
		"""
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

#take sampleName as input
#return path to peak examples and heatmap png
def get_peak_example(sampleName):
	peakMode = get_peak_mode(sampleName)
	if peakMode == "histone":
		return sampleDir + "/" + sampleName + "/HomerPeak.histone/peak.examples.png"
	elif peakMode == "factor":
		return sampleDir + "/" + sampleName + "/HomerPeak.factor/peak.examples.png"

def get_heatmap(sampleName):
	peakMode = get_peak_mode(sampleName)
	if peakMode == "histone":
		return sampleDir + "/" + sampleName + "/HomerPeak.histone/heatmap.exBL.png"
	elif peakMode == "factor":
		return sampleDir + "/" + sampleName + "/HomerPeak.factor/heatmap.exBL.1rpm.png"

rule create_report_per_sample:
	input:
		alnStat = qcDir + "/alignStat.txt",
		uniqFragCnt = qcDir + "/uniqFragCnt.txt",
		baseFreqPNG = sampleDir + "/{sampleName}/QC/base_freq.png",
		fragDistPNG = sampleDir + "/{sampleName}/QC/fragLen.dist.png",
		fragQC = sampleDir + "/{sampleName}/QC/fragMix.txt",
		peakExample = lambda wildcards: get_peak_example(wildcards.sampleName),
		heatmap = lambda wildcards: get_heatmap(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/QC/Report.html"
	message:
		"Creating report for sample... [{wildcards.sampleName}]"
	params:
		group = lambda wildcards: get_group(wildcards.sampleName)
	shell:
		"""
		module load Cutlery/1.0
		module load ImageMagick/6.9.12
		cnr.createSampleReportHTML.r -o Report -g {params.group} -s {sampleDir}/{wildcards.sampleName} -q {qcDir} -f fragMix
		"""


rule create_final_report:
	input:
		uniqFragCnt = qcDir + "/uniqFragCnt.txt",
		histPeakExamples = expand(sampleDir + "/{sampleName}/HomerPeak.histone/peak.examples.png", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		tfPeakExamples = expand(sampleDir + "/{sampleName}/HomerPeak.factor/peak.examples.png", sampleName = samples.Name[samples.PeakMode=="factor"].tolist()),
		fragDist = expand(sampleDir + "/{sampleName}/QC/fragLen.dist.txt", sampleName=samples.Name.tolist()),
		fragQC = expand(sampleDir + "/{sampleName}/QC/fragMix.txt", sampleName=samples.Name.tolist()),
		histoneHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.histone/heatmap.exBL.png", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		factorHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.factor/heatmap.exBL.1rpm.png", sampleName = samples.Name[samples.PeakMode=="factor"].tolist())
	output:
		"Report.html"
	message:
		"Creating final report in HTML..."
	shell:
		"""
		module load Cutlery/1.0
		module load ImageMagick/6.9.12
		cnr.createReportHTML.r -o Report -s {sampleDir} -q {qcDir} -f fragMix
		"""

rule create_report_per_sample_pooled:
	input:
		uniqFragCnt = qcDir + "/uniqFragCnt.txt",
		baseFreqPNG = sampleDir + "/{sampleName}/QC/base_freq.png",
		fragDistPNG = sampleDir + "/{sampleName}/QC/fragLen.dist.png",
		fragQC = sampleDir + "/{sampleName}/QC/fragMix.txt",
		peakExample = lambda wildcards: get_peak_example(wildcards.sampleName),
		heatmap = lambda wildcards: get_heatmap(wildcards.sampleName)
	output:
		sampleDir + "/{sampleName}/QC/Report_pooled.html"
	message:
		"Creating report for sample... [{wildcards.sampleName}]"
	params:
		group = lambda wildcards: get_group(wildcards.sampleName)
	shell:
		"""
		module load Cutlery/1.0
		module load ImageMagick/6.9.12
		cnr.createPooledSampleReportHTML.r -o Report_pooled -g {params.group} -s {sampleDir}/{wildcards.sampleName} -q {qcDir}		
		"""

rule create_final_report_pooled:
	input:
		histPeakExamples = expand(sampleDir + "/{sampleName}/HomerPeak.histone/peak.examples.png", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		tfPeakExamples = expand(sampleDir + "/{sampleName}/HomerPeak.factor/peak.examples.png", sampleName = samples.Name[samples.PeakMode=="factor"].tolist()),
		fragDist = expand(sampleDir + "/{sampleName}/QC/fragLen.dist.txt", sampleName=samples.Name.tolist()),
		fragQC = expand(sampleDir + "/{sampleName}/QC/fragMix.txt", sampleName=samples.Name.tolist()),
		histoneHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.histone/heatmap.exBL.png", sampleName = samples.Name[samples.PeakMode=="histone"].tolist()),
		factorHeatmap = expand(sampleDir + "/{sampleName}/HomerPeak.factor/heatmap.exBL.1rpm.png", sampleName = samples.Name[samples.PeakMode=="factor"].tolist())
	output:
		"Report_pooled.html"
	message:
		"Creating final report in HTML..."
	shell:
		"""
		module load Cutlery/1.0
		module load ImageMagick/6.9.12
		cnr.createPooledReportHTML.r -o Report_pooled -s {sampleDir} -q {qcDir}
		"""
