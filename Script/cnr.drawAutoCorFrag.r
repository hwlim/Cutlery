#!/usr/bin/env Rscript


#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
#suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
suppressPackageStartupMessages(library('tools', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))

source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))


# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix including path, default=<same with the src file excluding an extension under current directorys>")
#	make_option(c("-d","--maxDist"), default=1000, help="Max fragment length, x-axis for plotting. default=1000")
#	make_option(c("-f","--bamFlag"), default="0x2", help="flag for bam records. NULL is allowed to unset. Ignored for bed file. default=0x2 (concordant pairs only)"),
#	make_option(c("-F","--bamUnFlag"), default="0x400", help="flag for bam records to exclude, NULL is allowed to unset. Ignored for bed file. default=0x400 (exclude duplicates)")
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <*.bed.gz or *.bed>",
	description="Description:
	Calculate fragment auto-correlation of a given fragment bed file,
	Draw a auto-correlation plot,
	Calculate A0-index and A1-index for fragment aggregation. A0 includes 0 distance / A1 exclude 0 distance.
Output:
	- <outPrefix>.acor.txt
	- <outPrefix>.acor.<pdf/png>",
	 option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires input file")
} else {
	src <- arguments$args[1]
}


# Option handling
opt=arguments$options

outPrefix=opt$outPrefix
#maxDist=opt$maxDist
maxDist=1999
assertFileExist(src)


if(FALSE){
	src="JP1010_CnR_Hnf4a.frag.bed.gz"
	maxDist=1000
	outPrefix=NULL	
}

if(is.null(outPrefix)){
	outPrefix = sub(".bed.gz|.bed","", src)
	outPrefix = basename(outPrefix)
}
desDir=dirname(outPrefix)


write(sprintf("Checking fragment auto-correlation"), stderr())


system(sprintf("mkdir -p %s", desDir))
des.acor = sprintf("%s.acor.txt", outPrefix)
des.pdf = sprintf("%s.acor.pdf", outPrefix)
des.png = sprintf("%s.acor.png", outPrefix)
	
write(sprintf("  - %s", src), stderr())
#cmd=sprintf("cnr.autoCorFrag.sh -o %s -m %d %s", des.acor, maxDist, src)
cmd=sprintf("cnr.autoCorFrag.sh -o %s -m homer %s", des.acor, src)
system(cmd)

data.acor = read.delim(des.acor, header=TRUE)
a0_index=sum(data.acor[1:201,2])/sum(data.acor[,2])*100
a1_index=sum(data.acor[2:201,2])/sum(data.acor[2:(maxDist+1),2])*100

g0 = ggplot(data=data.acor, aes(x=Distance, y=NormFreq)) + geom_line() +
	labs(title=sprintf("Fragment Auto-Correlation: A0-index = %.01f", a0_index),
		x="Distance to the Closest Next Fragment (bp)",
		y="Normalized Frequency") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g1 = ggplot(data=data.acor[2:nrow(data.acor),], aes(x=Distance, y=NormFreq)) + geom_line() +
	labs(title=sprintf("Fragment Auto-Correlation excl. 0: A1-index = %.01f", a1_index),
		x="Distance to the Closest Next Fragment (bp)",
		y="Normalized Frequency") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g0.log = ggplot(data=data.acor, aes(x=Distance, y=log2(NormFreq))) + geom_line() +
	labs(title=sprintf("Fragment Auto-Correlation: A0-index = %.01f", a0_index),
		x="Distance to the Closest Next Fragment (bp)",
		y="log2(Normalized Frequency)") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g1.log = ggplot(data=data.acor[2:nrow(data.acor),], aes(x=Distance, y=log2(NormFreq))) + geom_line() +
	labs(title=sprintf("Fragment Auto-Correlation exc. 0: A1-index = %.01f", a1_index),
		x="Distance to the Closest Next Fragment (bp)",
		y="log2(Normalized Frequency)") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g.linear = plot_grid(g0, g1, nrow=2)
g.log = plot_grid(g0.log, g1.log, nrow=2)
g = plot_grid(g.linear, g.log, nrow=1)
ggsave(des.pdf, g, width=10, height=5)
ggsave(des.png, g, width=10, height=5, dpi=200)

