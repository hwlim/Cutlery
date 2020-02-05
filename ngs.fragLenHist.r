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
	make_option(c("-m","--maxLen"), default=1000, help="Max fragment length, x-axis for plotting. default=1000")
#	make_option(c("-f","--bamFlag"), default="0x2", help="flag for bam records. NULL is allowed to unset. Ignored for bed file. default=0x2 (concordant pairs only)"),
#	make_option(c("-F","--bamUnFlag"), default="0x400", help="flag for bam records to exclude, NULL is allowed to unset. Ignored for bed file. default=0x400 (exclude duplicates)")
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <bam or bed.gz>",
	description="Description:
	Check and visualize fragment length distribution for a paired-end BAM file. Considers chromosomes starting with \"chr\" only
Input:
	Paired-end BAM file or fragment bed file
Output:
	- <outPrefix>.dist.txt
	- <outPrefix>.dist.png",
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

#xLabel=opt$xlabel
#yLabel=opt$ylabel
#mainTitle=opt$title

outPrefix=opt$outPrefix
#bamFlag=opt$bamFlag
#bamUnFlag=opt$bamUnFlag
maxLen=opt$maxLen
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
	write(sprintf("Input is BED format"), stderr())
	mode="bed"
}else if(ext=="bam"){
	write(sprintf("Input is BAM format"), stderr())
	mode="bam"
}else{
	stop(sprintf("Invalid file type: %s", ext))
}

#outPrefix=opt$outPrefix
#if( is.null(outPrefix) ) outPrefix=basename

write(sprintf("Checking fragment length distribution"), stderr())
#if( mode=="bam" ){
#	if(bamFlag=="NULL"){
#		bamFlag=""
#	}else{
#		bamFlag=sprintf("-f %s", bamFlag)
#		write(sprintf("  Filtering by SAM flag: %s", bamFlag), stderr())
#	}
#
#	if(bamUnFlag=="NULL"){
#		bamUnFlag=""
#	}else{
#		bamUnFlag=sprintf("-F %s", bamUnFlag)
#		write(sprintf("  Filtering by SAM flag: %s", bamUnFlag), stderr())
#	}
#	if( bamFlag=="" && bamUnFlag=="" ) write(sprintf("  No filtering by SAM flag"), stderr())
#	write("", stderr())
#}

system(sprintf("mkdir -p %s", desDir))
des.dist = sprintf("%s.dist.txt", outPrefix)
des.hist = sprintf("%s.dist.png", outPrefix)
	
write(sprintf("  - %s", src), stderr())
if( mode=="bam" ){
	cmd=sprintf("bamToBed -bedpe -i %s 2>&1 | grep ^chr | gawk 'BEGIN{printf \"fragLen\\tCnt\\n\"; maxLen=%d }{ if($1==\".\" || $3==\".\") next; d=$6-$2; if(d>maxLen){d=maxLen}; cnt[d]++ }END{ for( i=1;i<=maxLen;i=i+1 ) printf \"%%d\\t%%d\\n\", i, cnt[i] }' > %s", src, maxLen, des.dist)
}else{
	cmd=sprintf("zcat %s | grep ^chr | gawk 'BEGIN{printf \"fragLen\\tCnt\\n\"; maxLen=%d }{ d=$3-$2; if(d>maxLen){d=maxLen}; cnt[d]++ }END{ for( i=1;i<=maxLen;i=i+1 ) printf \"%%d\\t%%d\\n\", i, cnt[i] }' > %s", src, maxLen, des.dist)
}
system(cmd)
data.dist = read.delim(des.dist, header=TRUE)

g1 = ggplot(data=data.dist, aes(x=fragLen, y=Cnt/1000)) + geom_line() +
	labs(title="Fragment Length Distribution", x="Fragment Length", y="Frequency (x1000)")
g2 = ggplot(data=data.dist, aes(x=fragLen, y=log10(Cnt))) + geom_line() +
	labs(title="Fragment Length Distribution", x="Fragment Length", y="log10(Frequency)")
g = plot_grid(g1, g2, nrow=2)
ggsave(des.hist, g, width=6, height=6)

