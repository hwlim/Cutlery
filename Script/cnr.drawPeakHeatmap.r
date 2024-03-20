#!/usr/bin/env Rscript

# Draw a bigwig heatmap for CUT&Run data
# Written by Hee Woong Lim & Christopher Ahn


suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('tools', quiet=TRUE))

source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))

# command line option handling
option_list <- list(
	make_option(c("-t","--title"), default="Peak Heatmap", help="Main Title. default=Peak Heatmap"),
	make_option(c("-w","--width"), default=2000, help="Visualization window, default=2000; Leave as 2000 for TF samples and use 10000 instead for histone samples."),
	make_option(c("-b","--binSize"), default=20, help="Bin size for tiled-average, default=20"),
	make_option(c("-m","--margin"), default="0,0.5,2,0.5", help="Outside margin (oma), default=0,0.5,2,0.5"),
	make_option(c("-s","--size"), default="3,3", help="Comma-separated figure size, xSize,ySize, in inch, default=3,3"),
	make_option(c("-q","--maxColor"), default="red", help="Color for the maximum value to use for color function, default=red"),
	make_option(c("-p","--dpi"), default=200, help="DPI for pdf->png conversion, default=200"),
	make_option(c("-d","--dataDir"), default="NULL", help="Directory containing extracted data. default=<outPrefix>.data"),
	make_option(c("-o","--outPrefix"), default="PeakHeatmap", help="Output prefix. default=PeakHeatmap")
)
parser <- OptionParser(usage = "%prog [options] <bed1> <bed2> ... <bw1> <bw2> ...", option_list=option_list,
			description = "Draw a matrix of CUT&Run data heatmaps for one or more bed files using given bigwig files
Input:
	- BED: One or more bed file to draw heatmaps
	- BigWig: ordered list of two or four bed files
		Two files: NFR and NUC bigwig files
		Four files: NFR / NFR IgG / NUC / NUC IgG
		Heatmap columns are named accordingly
Output:
	- a matrix of heatmaps where bed bed file corresponds to a row, and each bigwig file corresponds to a column.
		Colorscale is automatically decided separately in NFR and NUC
		If IgG bigwig files are given, they are visualized with the same color scale matching NFR and NUC
			e.g. NFR & NFR IgG / NUC & NUC IgG
		If more than one bed file is given, row names are annotated for each bed following the file name excluding the extension
Note:
	- Not appropriate for 1bp-sensitive jobs.
	- If data files already exist, only the plots are re-drawn."
	)

arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	q()
} else {
	srcL=arguments$args
}

# Option handling
opt=arguments$options

mainTitle=opt$title
width=opt$width
binSize=opt$binSize
axisLim=opt$axisLim
size=opt$size
dpi=opt$dpi
margin=opt$margin
dataDir=opt$dataDir
outPrefix=opt$outPrefix
colName=opt$colName
maxColor = opt$maxColor

stopifnot(isColor(maxColor))

###################################
## Input file classification
srcL.bed=NULL
srcL.bw=NULL
for( src in srcL ){
	ext=file_ext(src)
	if(ext=="bed"){
		srcL.bed = c(srcL.bed, src)
	}else{
		if(ext=="bw"){
			assertFileExist(src)
		}
		srcL.bw = c(srcL.bw, src)
	}
}
assertFileExist(srcL.bed)

###################################
## Option processing & validation

N.bed = length(srcL.bed)
N.bw = length(srcL.bw)
stopifnot( N.bed > 0 && N.bw > 0 )

#Define column names and re-order bw files if ctrl sample exists
if( N.bw == 4 ){
	nameL.bw = c("NFR", "NFR_Ctrl", "NUC", "NUC_Ctrl")
	margin = "0,3.8,2.5,3.8"
}else if ( N.bw == 2 ) {
    nameL.bw = c("NFR", "NUC")
	margin = "0,3.8,1.5,3.8"
} else {
	stop("Number of bigwig file should either be 2 or 4")
}

nameL.bed = sapply(srcL.bed, function(x) sub(".bed$","",basename(x)))

width=as.numeric(width)
stopifnot(is.numeric(width))

binSize=as.numeric(binSize)
stopifnot(is.numeric(binSize))

size=as.numeric(strsplit(size,",")[[1]])
stopifnot(length(size)==2)
stopifnot(is.numeric(size))

margin=as.numeric(strsplit(margin,",")[[1]])
stopifnot(length(margin)==4)
stopifnot(is.numeric(margin))

stopifnot(is.numeric(dpi) && dpi > 0)

if( dataDir=="NULL" ){
	dataDir=sprintf("%s.data", outPrefix)
}
system(sprintf("mkdir -p %s", dataDir))

####################################
## BigWig data extraction & read
dataL=list()
for( i in 1:N.bw ){
	bwName=nameL.bw[i]
	dataL[[bwName]]=list()
	bw=srcL.bw[i]
	
	for( j in 1:N.bed){
		bedName=nameL.bed[j]
		bed=srcL.bed[j]
		tmp.bed = read.delim(bed, header=FALSE, stringsAsFactors=FALSE)

		if(file.exists(bw)){
			## unstranded / single bw file
			des=sprintf("%s/%s.%s.%dbp.gz", dataDir, bedName, bwName, width)
			if(file.exists(des)){
				write(sprintf("Warning: %s already exists, pass", des), stderr())
				data = read.delim(gzfile(des), header=FALSE, row.names=1)
			}else{
				write(sprintf("Extracting => %s", des), stderr())
				write(sprintf("  - Bed: %s", bed), stderr())
				write(sprintf("  - BigWig: %s", bw), stderr())
				data = extractBigWigData(bed, bw, width=width, binSize=binSize, des)
			}
		}
		dataL[[bwName]][[bedName]]=data
	}
}

## Color scheme maximum
axisLim=NULL
for( bwName in names(dataL) ){
    tmp=do.call(rbind, dataL[[bwName]])
    maxVal=median(apply(abs(tmp), 1, max))
    if(maxVal==0) mean(apply(abs(tmp), 1, max))
    if(maxVal==0) maxVal=0.1
    axisLim=c(axisLim, maxVal)
}

if (N.bw == 4){
    axisLim[2] = axisLim[1]
    axisLim[4] = axisLim[3]
}

colorFun=function(n){
	return(colorpanel(n, "blue","white", maxColor))
}

pdf(sprintf("%s.pdf", outPrefix), width=size[1], height=size[2], bg="white")
par(oma=margin)
drawHeatmapList2(dataL, maxVal=axisLim, colorFun=colorFun, margin=c(2,0.2,3,0.2), cex=0.5)
mtext(paste0(mainTitle, "\n", width, "bp"), outer=TRUE, cex=0.7)
dev.off()
system(sprintf("convert -density %d %s.pdf %s.png", dpi, outPrefix, outPrefix))
system(sprintf("rm %s.pdf", outPrefix))