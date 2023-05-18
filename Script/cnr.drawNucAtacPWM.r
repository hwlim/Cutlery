#!/usr/bin/env Rscript


#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
#suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
suppressPackageStartupMessages(library('tools', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('reshape2', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))

source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))


# command line option handling
option_list <- list(
	make_option(c("-o","--out"), default=NULL, help="Input filen name + pdf")
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
des=opt$out


if(FALSE){
	src="LBP_FOXA2.CnR.PWM.txt"
	out="LBP_FOXA2.CnR.PWM.pdf"
}


if(is.null(des)) des = sprintf("%s.pdf", sub(".txt$","", basename(src)))

data=read.delim(src,header=FALSE, skip = 9)
rownames(data) = c("A","C","G","T")
colnames(data) = 1:ncol(data)
data = t(data)

df = melt(data)
colnames(df)=c("Pos","Base","Value")
g = ggplot(df, aes(x=Pos, y=Value, colour=Base)) + geom_line() + geom_point() + ggtitle(src)
ggsave(des, g, width=6, height=3)