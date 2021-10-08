#!/usr/bin/env Rscript


#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
#suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
suppressPackageStartupMessages(library('reshape', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))

source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))


# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix including path. Required"),
	make_option(c("-m","--maxDist"), default=1000, help="Max fragment length, x-axis for plotting. default=1000")
#	make_option(c("-f","--bamFlag"), default="0x2", help="flag for bam records. NULL is allowed to unset. Ignored for bed file. default=0x2 (concordant pairs only)"),
#	make_option(c("-F","--bamUnFlag"), default="0x400", help="flag for bam records to exclude, NULL is allowed to unset. Ignored for bed file. default=0x400 (exclude duplicates)")
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <a list of fragment files> or <a file containing fragment file list & annotation>",
	description="Description:
	For a given list of fragment bed files:
	- Calculate fragment auto-correlation
	- Draw all-in-one auto-correlation plot sorted by A0- or A1-index
	- Print the results to file
Input:
	1) a list of more than one fragment files
	2) a file containing fragment file list and annotation
		Possibly can have three columns: Fragment / Group / Name
		Fragment: input fragment files, required
		Group: Optional. If present, auto-correlation plot is color-coded by group
			If not, they are sorted by A0/A1-index and color-coded sequentially
		Name: Optional. alias names to print for input files
			If not, fragment files names are used instead
Output:
	- <outPrefix>.acor.txt: columns are auto-correlation for input files
	- <outPrefix>.aindex.txt: input annotation and A0-/A1-index columns
	- <outPrefix>.acor.<pdf/png>",
	 option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires input file")
} else {
	srcL <- arguments$args
}


# Option handling
opt=arguments$options
outPrefix=opt$outPrefix
maxDist=opt$maxDist
assertFileExist(src)


if(FALSE){
	srcL="JP1010_CnR_Hnf4a.frag.bed.gz"
	maxDist=1000
	outPrefix=NULL	
}

if(length(srcL) == 1){
	data.src = read.delim(srcL, header=TRUE, stringsAsFactors=FALSE)
	stopifnot("Fragment" %in% colnames(data.src))
	srcL=data.src$Fragment
	if("Group" %in% colnames(data.src)){
		groupL = data.src$Group
	}else{
		groupL=NULL
	}
	if("Name" %in% colnames(data.src)){
		nameL = data.src$Name
	}else{
		nameL = NULL
	}
}else{
	groupL=NULL
	nameL = srcL
}
assertFileExit(srcL)


if(is.null(outPrefix)) stop("outPrefix (-o) must be specified")
desDir=dirname(outPrefix)







system(sprintf("mkdir -p %s", desDir))
des.acor = sprintf("%s.acor.txt", outPrefix)
des.aindex = sprintf("%s.aindex.txt", outPrefix)
des.pdf = sprintf("%s.acor.pdf", outPrefix)
des.png = sprintf("%s.acor.png", outPrefix)

des.tmp = tempfile()
data.acor=NULL
data.aindex=NULL
write(sprintf("Checking fragment auto-correlation of multiple samples"), stderr())
for(src in srcL){
	write(sprintf("  - Processing %s", src), stderr())
	cmd=sprintf("cnr.autoCorFrag.sh -o %s -m %d %s", des.tmp, maxDist, src)
	system(cmd)

	tmp = read.delim(des.acor, header=TRUE)
	a0_index=sum(tmp[1:201,2])/sum(tmp[,2])*100
	a1_index=sum(tmp[2:201,2])/sum(tmp[2:(maxDist+1),2])*100

	if(is.null(data.acor)){
		data.acor = tmp
		data.aindex = c(a0_index, a1_index)
	}else{
		data.acor = cbind(data.acor, tmp[,2])
		data.aindex = cbind(data.aindex, c(a0_index, a1_index))
	}
}

colnames(data.acor) = c("Distance", nameL)
data.aindex = t(data.aindex)
colnames(data.aindex) = c("A0_Index","A1_Index")
data.aindex = data.frame(Name = nameL, data.aindex)
if(!is.null(groupL)) data.aindex = data.frame(data.aindex, Group = groupL)

write.table(data.acor, des.acor, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(data.aindex, des.aindex, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

acor.melt = melt(data.acor, id.var="Distance")

g0 = ggplot(data=data.acor, aes(x=Distance, y=value, group=variable)) + geom_line() +
	labs(title="Fragment Auto-Correlation",
		x="Distance to the Closest Next Fragment (bp)",
		y="Normalized Frequency") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g1 = ggplot(data=data.acor[2:nrow(data.acor),], aes(x=Distance, y=NormFreq)) + geom_line() +
	labs(title="Fragment Auto-Correlation excl. 0",
		x="Distance to the Closest Next Fragment (bp)",
		y="Normalized Frequency") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g0.log = ggplot(data=data.acor, aes(x=Distance, y=log2(NormFreq))) + geom_line() +
	labs(title="Fragment Auto-Correlation",
		x="Distance to the Closest Next Fragment (bp)",
		y="log2(Normalized Frequency)") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g1.log = ggplot(data=data.acor[2:nrow(data.acor),], aes(x=Distance, y=log2(NormFreq))) + geom_line() +
	labs(title="Fragment Auto-Correlation exc. 0",
		x="Distance to the Closest Next Fragment (bp)",
		y="log2(Normalized Frequency)") +
	geom_vline(xintercept=200, color="purple", linetype="dashed")

g.linear = plot_grid(g0, g1, nrow=2)
g.log = plot_grid(g0.log, g1.log, nrow=2)
g = plot_grid(g.linear, g.log, nrow=1)
ggsave(des.pdf, g, width=10, height=5)
ggsave(des.png, g, width=10, height=5, dpi=200)

