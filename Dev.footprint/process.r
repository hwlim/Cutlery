#!/usr/bin/env Rscript

#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
#suppressPackageStartupMessages(library('robustbase', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/motifR.r", Sys.getenv("COMMON_LIB_BASE")))

# command line option handling
option_list <- list(
	make_option(c("-x","--xlabel"), default="xLabel", help="X-axis label"),
	make_option(c("-y","--ylabel"), default="yLabel", help="Y-axis label"),
	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
	make_option(c("-o","--outfile"), default="barplot.png", help="Output file, with .png extension. [default: barplot.png]")
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] dataFile", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires a data file")
} else {
	dataFile <- arguments$args[1]
}

# Option handling
opt=arguments$options

xLabel=opt$xlabel
yLabel=opt$ylabel
mainTitle=opt$title


if(FALSE){
	## Homer motif dir
	src.motifDir="TestData/hESC_Sox2-HomerPeak.factor-peak.homer.exBL.1rpm.bed.all.noBG/"
	src.peak="TestData/hESC_Sox2-HomerPeak.factor-peak.homer.exBL.1rpm.bed"
	minTargetPercent = 10
	maxMotifCount = 10
	genome="hg19"
	bwPrefix="TestData/hESC_Sox2"
}

motifDir=sprintf("%s/homerResults", motifDir)
assertDirExist(motifDir)
assertFileExist(src.peak)
assertFileExist(sprintf("%s.%s.bw", bwPrefix, c("plus","minus")))

## Parse homer motif header line (1st line) and return following information
## - Best guess name
## - Motif score threshold
## - log(p-value)
## - p-value
## - Target %
## - Background %
parseHomerMotifHeader=function(hdr){
	element = strsplit(hdr, "\t")[[1]]
	result=list()

	## Usual format: 1-AWRACAAWRG,BestGuess:SOX15/MA1152.1/Jaspar(0.974)
	## Special formats:
	##		2-TTGTTATGCAAA,BestGuess:Pou5f1::Sox2/MA0142.1/Jaspar(0.933)
	##		4-NRCATTCCTT,BestGuess:TEAD(TEA)/Fibroblast-PU.1-ChIP-Seq(Unpublished)/Homer(0.962)
	##		==> :: vs :
	##		==> () in the name
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
## Motif selection
src.motifCollection="./selectedMotif.motif"
append=FALSE
motifNameSelectL=NULL
for( i in 1:maxMotifCount ){
	# i=2
	src.motif = sprintf("%s/motif%d.motif", motifDir, i)
	if( ! file.exists(src.motif) ) break
	write(sprintf("Processing motif #%d\n  - %s", i, src.motif), stderr())

	hdr = readLines(src.motif, n=1)
	hdr.parsed = parseHomerMotifHeader(hdr)
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

###############################
## Motif scan

## Retain motif loci only; Discard peak annotation result
src.motifLoci = "./motifScan.bed"
cmd=sprintf("annotatePeaks.pl %s %s -m %s -noann -nogene -mbed %s > /dev/null", src.peak, genome, src.motifCollection, src.motifLoci)
system(cmd)

data.scan = read.delim(src.motifLoci, header=FALSE, skip=1, stringsAsFactors=FALSE)
motifL=unique(data.scan[,4])
N.motif = length(motifL)




####################################
## Filtering by footprint contrast

## Function to measure footprint contrast
## i.e. MNase digestion within anchor region vs flanking margins
getFootprintContrast=function(anchor, margin){

	anchor.uniqName = data.frame( anchor[,1:3], sprintf("anchor.%d", 1:nrow(anchor)), anchor[,c(5,6,4)] )
	anchor.ext = anchor.uniqName
	anchor.ext[,2] = bed.anchor[,2] - margin
	anchor.ext[,3] = bed.anchor[,3] + margin
	anchor.ext[,5] = 0
	anchor.ext[,6] = "+"


	src.tmpBed = tempfile()
	write.table(anchor.ext[,1:6], src.tmpBed, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
	write(sprintf("\tplus = %s", des.bwPlus), stderr())
	cmd = sprintf("bwtool extract bed -fill=0 -decimals=5 -tabs %s %s /dev/stdout | cut -f 4,8- | gzip -c > %s", src.tmpBed, src.bwPlus, des.bwPlus)
	system(cmd)
	write(sprintf("\tminus = %s", des.bwMinus), stderr())
	cmd = sprintf("bwtool extract bed -fill=0 -decimals=5 -tabs %s %s /dev/stdout | cut -f 4,8- | gzip -c > %s", src.tmpBed, src.bwMinus, des.bwMinus)
	system(cmd)
	unlink(src.tmpBed)
	profileL = readExoProfile2(srcPlus=des.bwPlus, srcMinus=des.bwMinus, direc=anchor.ext[,6])


	idx.anchor = (margin + 1):(margin+anchorSize)
	idx.flank = c( 1:margin, (margin+anchorSize+1):(2*margin+anchorSize) )

	signal.anchor = apply( profileL[["plus"]][,idx.anchor], 1, mean ) + abs(apply( profileL[["minus"]][,idx.anchor], 1, mean ))
	signal.flank = apply( profileL[["plus"]][,idx.flank], 1, mean ) + abs(apply( profileL[["minus"]][,idx.flank], 1, mean ))


	### All data
	colnames(anchor) = c("Chr","Start","End","Name","Null","Direc")
	data.final = data.frame(anchor,
					Anchor = signal.anchor,
					Flank = signal.flank,
					Contrast = log2( (signal.flank + pseudo) / (signal.anchor + pseudo) )
				)
}


margin=20
for( i in 1:N.motif ){
	motifName =  motifL[i]
	src.anchor = sprintf("%s.scan.bed", motifName)

	anchor = data.scan[data.scan[,4]==motifName,]
	write.table(bed.select, src.anchor, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

	## Contrast
	contrast = getFootprintContrast(anchor, margin)
	write.table(contrast, src.anchor, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

	##
}


####################################
## Footprint visualization
