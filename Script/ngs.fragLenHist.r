#!/usr/bin/env Rscript


#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
#suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
suppressPackageStartupMessages(library('tools', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('plotly', quiet=TRUE))
require(reshape2)

source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))


# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix including path, default=<same with the src file excluding an extension under current directorys>"),
	make_option(c("-l","--maxLen"), default=1000, help="Max fragment length, x-axis for plotting. default=1000"),
	make_option(c("-n","--name"), default=NULL, help="Sample name to display at top. default=<input file name>"),
	make_option(c("-i","--interactive"), default=FALSE, action="store_true", help="If set, interactive plotly plot is also generated in html"),
	make_option(c("-c","--cSorted"), default=FALSE, action="store_true", help="If set, input bam file is assumed to be coordinate-sorted not name-sorted"),
	make_option(c("-a","--getAverage"), default=FALSE, action="store_true", help="IF set, return average fragment length by WMA will be returned to STDOUT and highlighted in the plot")
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <bam or bed.gz>",
	description="Description:
	Check and visualize fragment length distribution for a paired-end BAM file or a fragment bed file.
	Considers chromosomes starting with \"chr\" only
Input:
	Paired-end BAM file or fragment bed file
Output:
	- <outPrefix>.txt
	- <outPrefix>.png
	- <outPrefix>.html (if -i is set)",
	 option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires bam file(s)")
} else {
	src <- arguments$args[1]
}


# Option handling
opt=arguments$options


outPrefix=opt$outPrefix
maxLen=opt$maxLen
name=opt$name
drawPlotly = opt$interactive
cSorted=opt$cSorted
getAverage=opt$getAverage
assertFileExist(src)


if(FALSE){
	src="1.1.Align/DE_H3K9me3_Control_A6/align.sortByName.bam"
	maxLen=1000
	outPrefix="test"
	cSorted=FALSE
	getAverage=TRUE
	name="test"
	drawPlotly=FALSE
}

if(is.null(name)) name=src

if(is.null(outPrefix)){
	outPrefix = sub(".bed.gz|.bam","", src)
	outPrefix = basename(outPrefix)
}
desDir=dirname(outPrefix)
ext=file_ext(src)

if(ext=="gz"){
	write(sprintf("Input is BED format"), stderr())
	mode="bed"
}else if(ext=="bam"){
	write(sprintf("Input is BAM format"), stderr())
	mode="bam"
}else{
	stop(sprintf("Invalid file type: %s", ext))
}


write(sprintf("Checking fragment length distribution"), stderr())


system(sprintf("mkdir -p %s", desDir))
des.dist = sprintf("%s.txt", outPrefix)
des.hist = sprintf("%s.png", outPrefix)
	
write(sprintf("  - %s", src), stderr())

## Checking fragment length distribution
if(cSorted){
	cmd=sprintf("ngs.getFragLenHist.sh -o %s -l %d -c %s", des.dist, maxLen, src)
}else{
	cmd=sprintf("ngs.getFragLenHist.sh -o %s -l %d %s", des.dist, maxLen, src)
}
system(cmd)
data.dist = read.delim(des.dist, header=TRUE)

## Caculate average fragment length by WMA
## i.e. finding fragment length with maximum frequency after WMA
data.wma = data.dist
data.wma[nrow(data.wma),2] = 0
data.wma[,2] = wma(data.wma[,2],5)
avg = round(which.max(data.wma[,2]))



## Draw density plots: linear and log scale
if( getAverage ){
	# print to STDOUT
	write(avg, stdout())
	# Plot
	data.plot = data.frame(data.dist, WMA=data.wma[,2])
	colnames(data.plot) = c("FragLen","Raw","WMA")
	data.melt = melt(data.plot, measure.vars=c("Raw","WMA"))
	colnames(data.melt) = c("FragLen","Category","Count")
	g1 = ggplot(data=data.melt, aes(x=FragLen, y=Count/1000, colour=Category)) +
		geom_line() +
		labs(title=sprintf("%s (avg = %d bp)", name, avg), x="Fragment Length", y="Frequency (x1000)") +
		geom_vline(xintercept=avg, linetype="dashed", color="purple", linewidth=0.5)

	g2 = ggplot(data=data.melt, aes(x=FragLen, y=log10(Count), colour=Category)) +
		geom_line() +
		labs(title=sprintf("%s (avg = %d bp)", name, avg), x="Fragment Length", y="Frequency (x1000)") +
		geom_vline(xintercept=avg, linetype="dashed", color="purple", linewidth=0.5)
}else{
	g1 = ggplot(data=data.dist, aes(x=fragLen, y=Cnt/1000)) +
		geom_line() +
		labs(title=name, x="Fragment Length", y="Frequency (x1000)")
	g2 = ggplot(data=data.dist, aes(x=fragLen, y=log10(Cnt))) +
		geom_line() +
		labs(title=name, x="Fragment Length", y="log10(Frequency)")
}
g = plot_grid(g1, g2, nrow=2)
ggsave(des.hist, g, width=6, height=6)






## Draw interactive density plots
if(drawPlotly){
	des.html = sprintf("%s.html", outPrefix)
	fig <- plot_ly(data.dist, x = ~fragLen, y = ~Cnt, type = 'scatter', mode = 'lines')
	fig = fig %>% layout(
						title = name,
						xaxis = list(title = "Fragment Length"),
						yaxis = list(title = "Frequency")
						)

	tmp=sprintf("%s.html", tempfile(tmpdir="."))
	htmlwidgets::saveWidget(fig, tmp)
	system(sprintf("mv %s %s", tmp, des.html))
}
