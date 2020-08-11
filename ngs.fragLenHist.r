#!/usr/bin/env Rscript


#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
#suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
suppressPackageStartupMessages(library('tools', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('plotly', quiet=TRUE))

source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))


# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix including path, default=<same with the src file excluding an extension under current directorys>"),
	make_option(c("-l","--maxLen"), default=1000, help="Max fragment length, x-axis for plotting. default=1000"),
	make_option(c("-n","--name"), default=NULL, help="Sample name to display at top. default=<input file name>"),
	make_option(c("-i","--interactive"), default=FALSE, action="store_true", help="If set, interactive plotly plot is also generated in html")
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
	- <outPrefix>.dist.png
	- <outPrefix>.dist.html (if -p is set)",
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
assertFileExist(src)


if(FALSE){
	src="Sox2.frag.bed.gz"
	maxLen=1000
	outPrefix=NULL	
	des.dist="fragLen.dist.txt"
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
des.dist = sprintf("%s.dist.txt", outPrefix)
des.hist = sprintf("%s.dist.png", outPrefix)
	
write(sprintf("  - %s", src), stderr())
#if( mode=="bam" ){
#	cmd=sprintf("bamToBed -bedpe -i %s 2>&1 | grep ^chr | gawk 'BEGIN{printf \"fragLen\\tCnt\\n\"; maxLen=%d }{ if($1==\".\" || $3==\".\") next; d=$6-$2; if(d>maxLen){d=maxLen}; cnt[d]++ }END{ for( i=1;i<=maxLen;i=i+1 ) printf \"%%d\\t%%d\\n\", i, cnt[i] }' > %s", src, maxLen, des.dist)
#}else{
#	cmd=sprintf("zcat %s | grep ^chr | gawk 'BEGIN{printf \"fragLen\\tCnt\\n\"; maxLen=%d }{ d=$3-$2; if(d>maxLen){d=maxLen}; cnt[d]++ }END{ for( i=1;i<=maxLen;i=i+1 ) printf \"%%d\\t%%d\\n\", i, cnt[i] }' > %s", src, maxLen, des.dist)
#}
cmd=sprintf("ngs.getFragLenHist.sh -o %s -l %d %s", des.dist, maxLen, src)
system(cmd)
data.dist = read.delim(des.dist, header=TRUE)

g1 = ggplot(data=data.dist, aes(x=fragLen, y=Cnt/1000)) + geom_line() +
	labs(title=name, x="Fragment Length", y="Frequency (x1000)")
g2 = ggplot(data=data.dist, aes(x=fragLen, y=log10(Cnt))) + geom_line() +
	labs(title=name, x="Fragment Length", y="log10(Frequency)")
g = plot_grid(g1, g2, nrow=2)
ggsave(des.hist, g, width=6, height=6)


if(drawPlotly){
	des.html = sprintf("%s.dist.html", outPrefix)
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