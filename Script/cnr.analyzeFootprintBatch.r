#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/motifR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="cnr_footprint", help="Output prefix, default=cnr_footprint"),
	make_option(c("-n","--name"), default="NULL", help="Data name, default=<peak file name without suffix>"),
	make_option(c("-t","--minTargetPercent"), default=10, help="Minimum target % for motif selection, default=10"),
	make_option(c("-m","--maxMotifCount"), default=10, help="Maximum motif count to consider, default=10"),
	make_option(c("--genome"), default=NULL, help="Genome for Homer. e.g. hg38 or mm10. If not specified, parameters from Homer motif directory (motifFindingParameters.txt) is retrieved."),
	make_option(c("-b","--bwPrefix"), default=NULL, help="BigWig file prefix assuming <prefix>.{plus,minus}.bw, Required"),
	make_option(c("-p","--parallel"), default=1, help="Number of thread to use for parallel processing for visualization, default=1")
)
parser <- OptionParser(usage = "%prog [options] <bed> <Homer motif dir>", option_list=option_list,
	description="Description:
	Perform CUT&RUN footprint analysis using given peak file and Homer motif result directory
Input:
	- BED file containing motif scan regions, not necessarily fixed size
	- Homer motif results directory
Output:
	- <outPrefix>.0.selectedMotif.motif	: selected motif for scan
	- <outPrefix>.1.motifScan.bed		: motif scan result
	- <outPrefix>.2.<motif name>.select.bed	: filtered motif scan by footprint contrast
	- <outPrefix>.3.<motif name>/CnR.<suffix>	: footprint visualization for the filtered regions including
		*.sorted.fa
		*.logo.<pdf/png>
		*.viz.png
		*.viz.avg.<pdf/png>
		*.sorted.bed
		*.sorted.fa
	- <outPrefix>.4.complete		: flag file to show that analysis is successful and complete")
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires a data file")
} else {
	src.peak = arguments$args[1]
	src.motifDir = arguments$args[2]
}

# Option handling
opt=arguments$options

outPrefix=opt$outPrefix
name=opt$name
minTargetPercent=opt$minTargetPercent
maxMotifCount=opt$maxMotifCount
genome=opt$genome
bwPrefix=opt$bwPrefix
parallel=opt$parallel


if(FALSE){
	src.peak="TestData/hESC_Sox2-HomerPeak.factor-peak.homer.exBL.1rpm.bed"
	src.motifDir="TestData/hESC_Sox2-HomerPeak.factor-peak.homer.exBL.1rpm.bed.all.noBG/"
	outPrefix="TestOut/footprint"
	name="hESC_SOX2"
	minTargetPercent = 10
	maxMotifCount = 10
	genome=NULL
	bwPrefix="TestData/hESC_Sox2"
	parallel=2
}

motifDir=sprintf("%s/homerResults", src.motifDir)
assertDirExist(motifDir)
assertFileExist(sprintf("%s/homerResults.html", src.motifDir))
assertFileExist(src.peak)
assertFileExist(sprintf("%s.%s.bw", bwPrefix, c("plus","minus")))

## If genome is not specified, retrieve it from homer motif parameters
if(is.null(genome)){
	assertFileExist(sprintf("%s/motifFindingParameters.txt", src.motifDir))
	tmp = read.delim(sprintf("%s/motifFindingParameters.txt", src.motifDir), header=FALSE, sep=" ", stringsAsFactors=FALSE)
	genome = tmp[1,4]
}

## number of core to use
library("foreach")
library("doParallel")
N.core = detectCores()
if(parallel > N.core){
	write(sprintf("Warning: only %d cores available, reduce parallel from %d to %d", N.core, parallel, N.core), stderr())
	parallel = N.core
}

desDir=dirname(outPrefix)
system(sprintf("mkdir -p %s", desDir))


###########################################
## Output file names including intermediate results
src.motifCollection = sprintf("%s.0.selectedMotif.motif", outPrefix)
src.motifLoci = sprintf("%s.1.motifScan.bed", outPrefix)
## <outPrefix>.2.##.<motifName>.select.bed




## Parse homer motif NATIVE header line (1st line) and return following information
## Native --> direct output of Homer findMotifsGenome.pl
## - Best guess name
## - Motif score threshold
## - log(p-value)
## - p-value
## - Target %
## - Background %
parseHomerMotifNativeHeader=function(hdr){
	element = strsplit(hdr, "\t")[[1]]
	result=list()

	## Usual Homer motif name field (2nd field) format examples: 
	##	1-AWRACAAWRG,BestGuess:SOX15/MA1152.1/Jaspar(0.974)
	##		Normal, simply cut by : and /
	##	2-TTGTTATGCAAA,BestGuess:Pou5f1::Sox2/MA0142.1/Jaspar(0.933)
	##		Need to handle :: vs :
	##	4-NRCATTCCTT,BestGuess:TEAD(TEA)/Fibroblast-PU.1-ChIP-Seq(Unpublished)/Homer(0.962)
	##		Ned to handle () in the name
	tmp = strsplit(element[2],",BestGuess:", fixed=TRUE)[[1]][2]
	result[["name"]] = sub("::","_", strsplit(tmp, "[/()]")[[1]][1])
	result[["score"]] = element[3]
	result[["logp"]]= element[4]
	tmp = strsplit(element[6],"[%:(),]")[[1]]
	stopifnot(length(tmp)==12)
	stopifnot(tmp[1] == "T")
	stopifnot(tmp[6] == "B")
	stopifnot(tmp[11] == "P")

	result[["pvalue"]] = tmp[12]
	result[["target"]] = tmp[3]
	result[["background"]]  = tmp[8]

	return(result)
}

###############################
## Select motif sets for footprint visualization from Homer motif directory
write(sprintf("Collect motifs to scan:\n\tN = %d\n\tMin. target %% > %d %%", maxMotifCount, minTargetPercent), stderr())
append=FALSE
motifNameSelectL=NULL
for( i in 1:maxMotifCount ){
	# i=2
	src.motif = sprintf("%s/motif%d.motif", motifDir, i)
	if( ! file.exists(src.motif) ) break
	write(sprintf("Processing motif #%d\n  - %s", i, src.motif), stderr())

	hdr = readLines(src.motif, n=1)
	hdr.parsed = parseHomerMotifNativeHeader(hdr)
	write(sprintf("  best guess = %s", hdr.parsed[["name"]]), stderr())
	write(sprintf("  score = %s", hdr.parsed[["score"]]), stderr())
	write(sprintf("  p-value = %s", hdr.parsed[["pvalue"]]), stderr())
	write(sprintf("  log(p-value) = %s", hdr.parsed[["logp"]]), stderr())
	write(sprintf("  target %% = %s", hdr.parsed[["target"]]), stderr())
	write(sprintf("  background %% = %s", hdr.parsed[["background"]]), stderr())

	if( as.numeric(hdr.parsed[["pvalue"]]) > 1e-10 ){
		write(sprintf("  Warning: Skipping because p-value = %s > 1e-10", hdr.parsed[["pvalue"]]), stderr())
		next
	}
	if( as.numeric(hdr.parsed[["target"]]) < minTargetPercent ){
		write(sprintf("  Warning: Skipping because target %% = %s < %d", hdr.parsed[["target"]], minTargetPercent), stderr())
		next
	}
	motifNameSelectL = c(motifNameSelectL, hdr.parsed[["name"]])

	## Write to motif collection
	write(sprintf(">Motif%02d\tMotif%02d_%s\t%s\t%s\t0", i, i, hdr.parsed[["name"]], hdr.parsed[["score"]], hdr.parsed[["logp"]]),
		src.motifCollection, append=append)
	if(!append) append=TRUE
	system(sprintf("cat %s | grep -v \"^>\" >> %s", src.motif, src.motifCollection))
}
write(sprintf(""), stderr())



###############################
## Motif scan
## Retain motif loci only; Discard peak annotation result
write(sprintf("Scanning motifs"), stderr())
cmd=sprintf("annotatePeaks.pl %s %s -m %s -noann -nogene -mbed %s -cpu %d > /dev/null", src.peak, genome, src.motifCollection, src.motifLoci, parallel)
ret = system(cmd)
if(ret > 0) stop("Homer motif scan filed")
data.scan = read.delim(src.motifLoci, header=FALSE, skip=1, stringsAsFactors=FALSE)
motifL = unique(data.scan[,4])
N.motif = length(motifL)





## Function to measure footprint contrast
## i.e. MNase digestion within anchor region vs flanking margins
## NOTE: This routin assums all the anchor regions are the same size
## PLAN:
##	- allow variable sized anchor windows
getFootprintContrast=function(anchor, margin, bwPrefix, pseudo=0.1){

	## check if 4th column contains unique id. If not, reasign
	if( length(unique(anchor[,4])) != nrow(anchor) ){
		anchor = data.frame( anchor[,1:3], sprintf("%s:%d", anchor[,4], 1:nrow(anchor)), anchor[,c(5,6,4)] )
	}
	## anchor window size check if all the same
	anchorSize=anchor[1,3] - anchor[1,2]
	if(any(anchor[,3]-anchor[,2] != anchorSize)) stop("Window size are not the same")

	anchor.ext = anchor
	anchor.ext[,2] = anchor[,2] - margin
	anchor.ext[,3] = anchor[,3] + margin
	anchor.ext[,5] = 0
	anchor.ext[,6] = "+"
	src.bwPlus=sprintf("%s.plus.bw", bwPrefix)
	src.bwMinus=sprintf("%s.minus.bw", bwPrefix)

	## Signal exraction
	profileL = extractBigWigDataStranded1bp(anchor.ext[,1:6], bwPrefix, desPrefix=NULL)
	stopifnot(all(rownames(profileL[[1]])==anchor[,4]))
	stopifnot(all(rownames(profileL[[2]])==anchor[,4]))
	idx.anchor = (margin + 1):(margin+anchorSize)
	idx.flank = c( 1:margin, (margin+anchorSize+1):(2*margin+anchorSize) )
	signal.anchor = apply( profileL[["plus"]][,idx.anchor], 1, mean ) + abs(apply( profileL[["minus"]][,idx.anchor], 1, mean ))
	signal.flank = apply( profileL[["plus"]][,idx.flank], 1, mean ) + abs(apply( profileL[["minus"]][,idx.flank], 1, mean ))

	## All footprint contrast data
	colnames(anchor) = c("Chr","Start","End","Name","Null","Direc")
	result = data.frame(anchor,
					Anchor = signal.anchor,
					Flank = signal.flank,
					Contrast = log2( (signal.flank + pseudo) / (signal.anchor + pseudo) )
				)
	return(result)
}

####################################
## Filtering by footprint contrast
## PLAN:
##	- parallelize

write(sprintf("Checking footprint contrast"), stderr())
margin=20
srcL.selectedAnchor=list()
cntL.anchor = NULL
for( i in 1:N.motif ){
	# i=1
	motifName =  motifL[i]
	src.selectedAnchor = sprintf("%s.2.%s.select.bed", outPrefix, motifName)

	write(sprintf("  - Checking %s", motifName), stderr())
	anchor = data.scan[data.scan[,4]==motifName,]
	anchor[,4] = sprintf("%s:%d", anchor[,4], 1:nrow(anchor))

	## Contrast
	contrast = getFootprintContrast(anchor, margin, bwPrefix, pseudo=0.1)
	anchor.select = contrast[contrast$Contrast > 0, ]
	if(nrow(anchor.select) > 0){
		write.table(anchor.select, src.selectedAnchor, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
		srcL.selectedAnchor[[motifName]] = src.selectedAnchor
	}
	cntL.anchor = c(cntL.anchor, nrow(anchor.select))
}


####################################
## Footprint visualization
write(sprintf("Visualizing footprint"), stderr())
stopifnot(all(names(srcL.selectedAnchor) %in% motifL))

## Command line generation
indexL.select = which(cntL.anchor > 0) #match(names(srcL.selectedAnchor), motifL)
cmdL=NULL
for( i in indexL.select ){
	motifName = motifL[i]
	src.selectedAnchor = srcL.selectedAnchor[motifName]
	#mainTitle= sprintf("%s (N=%d)", motifName, cntL.anchor[i])
	vizPrefix=sprintf("%s.3.%s/CnR", outPrefix, motifName)
	cmd=sprintf("idom.visualizeExoBed.r -o %s -g %s -t %s -f -s %s %s", vizPrefix, genome, motifName, src.selectedAnchor, bwPrefix)
	cmdL = c(cmdL, cmd)
}

## EXecute commands
if(parallel > 1){
	## Parallel processing
	write(sprintf("Processing in parallel: N=%d", parallel), stderr())
	registerDoParallel(parallel)  # use multicore, set to the number of our cores
	ret = foreach( i=1:length(cmdL) ) %dopar% {
		write(sprintf("Visualizing %s", motifL[indexL.select[i]]), stderr())
		write(sprintf("\t%s", cmdL[i]), stderr())
		system(cmdL[i])
	}
	retL = unlist(ret)
}else{
	## Serial processing
	write(sprintf("Processing in serial"), stderr())
	retL = NULL
	for( i in 1:length(cmdL) ){
		write(sprintf("Visualizing %s", motifL[indexL.select[i]]), stderr())
		write(sprintf("\t%s", cmdL[i]), stderr())
		ret = system(cmdL[i])
		retL = c(retL, ret)
	}
}

## Write a flag file to tell analysis complete
if(any(retL > 0)){
	indL.fail = which(retL > 0)
	write("Error: Footprint visualization failed for:", stderr())
	for( i in indL.fail ) write(sprintf("  - %s", motifL[i]), stderr())
	q(1)
}else{
	write("Footprint visualization complete", stderr())
	write("", sprintf("%s.4.complete", outPrefix))
}
