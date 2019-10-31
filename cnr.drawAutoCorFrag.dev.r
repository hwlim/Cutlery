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
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix including path, default=<same with the src file excluding an extension under current directorys>"),
	make_option(c("-m","--maxDist"), default=1000, help="Max fragment length, x-axis for plotting. default=1000")
#	make_option(c("-f","--bamFlag"), default="0x2", help="flag for bam records. NULL is allowed to unset. Ignored for bed file. default=0x2 (concordant pairs only)"),
#	make_option(c("-F","--bamUnFlag"), default="0x400", help="flag for bam records to exclude, NULL is allowed to unset. Ignored for bed file. default=0x400 (exclude duplicates)")
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <*.bed.gz or *.bed>",
	description="Description:
	Calculate fragment auto-correlation of a given fragment bed file
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
maxDist=opt$maxDist
assertFileExist(src)


if(FALSE){
	src="Sox2.frag.bed.gz"
	maxLen=1000
	outPrefix=NULL	
}

if(is.null(outPrefix)){
	outPrefix = sub(".bed.gz|.bam","", src)
	outPrefix = basename(outPrefix)
}
desDir=dirname(outPrefix)
ext=file_ext(src)

if(ext=="gz"){
	write(sprintf("Input is gzipped bed format"), stderr())
	mode="bed.gz"
}else if(ext=="bed"){
	write(sprintf("Input is bed format"), stderr())
	mode="bed"
}else{
	stop(sprintf("Invalid file type: %s", ext))
}

#outPrefix=opt$outPrefix
#if( is.null(outPrefix) ) outPrefix=basename

write(sprintf("Checking fragment length distribution"), stderr())


system(sprintf("mkdir -p %s", desDir))
des.acor = sprintf("%s.acor.txt", outPrefix)
des.pdf = sprintf("%s.acor.pdf", outPrefix)
des.png = sprintf("%s.acor.png", outPrefix)
	
write(sprintf("  - %s", src), stderr())
cmd=sprintf("cnr.autoCorFrag.sh -o %s -m %d", des.acor, maxDist)
system(cmd)

data.acor = read.delim(des.acor, header=TRUE)
a_index0=sum(data.acor[1:201,2])/sum(data.acor[,2])*100
a_index1=sum(data.acor[2:201,2])/sum(data.acor[2:(maxDist+1),2])*100

g1 = ggplot(data=data.acor, aes(x=Distance, y=NormFreq)) + geom_line() +
	labs(title=sprintf("Fragment Auto-Correlation: A-index0 = %.01f / A-index1 = %.01f", a_index0, a_index1),
		x="Distance to the Closest Next Fragment (bp)",
		y="Normalized Frequency")
		
g2 = ggplot(data=data.dist, aes(x=fragLen, y=log10(Cnt))) + geom_line() +
	labs(title="Fragment Length Distribution", x="Fragment Length", y="log10(Frequency)")
g = plot_grid(g1, g2, nrow=2)
ggsave(des.hist, g, width=6, height=6)

