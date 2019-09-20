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
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix including path. required")
#	make_option(c("-f","--bamFlag"), default="0x2", help="flag for bam records. NULL is allowed to unset. Ignored for bed file. default=0x2 (concordant pairs only)"),
#	make_option(c("-F","--bamUnFlag"), default="0x400", help="flag for bam records to exclude, NULL is allowed to unset. Ignored for bed file. default=0x400 (exclude duplicates)")
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <spikeCnt.txt> <spikeCnt.txt> ...",
	description="Description:
	Collect target tag count and spikein count information from multiple files and make a table file and visualize bar plots
	Assuming each input file is in the format:
		Cartegory  Count
		Main       ####
		Spikein    ####
		Unknown    ####
Output:
	<outPrefix>.txt: spikein count table file with columns, Name / Main / Spikein
	<outPrefix>.png: barplots for read counts, spikein %, main-to-spikein ratio",
	 option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires input file(s)")
} else {
	srcL <- arguments$args
}


# Option handling
opt=arguments$options
outPrefix=opt$outPrefix

if( FALSE ){
	srcL=c(
		"DE_H3K9me3_Control_A6.spikeCnt.txt",
		"DE_H3K9me3_Control_H5.spikeCnt.txt",
		"DE_H3K9me3_FOXA123KD_A6.spikeCnt.txt",
		"DE_H3K9me3_FOXA123KD_H5.spikeCnt.txt"
	)
	outPrefix="spikein"
}
assertFileExist(srcL)
if(is.null(outPrefix)) stop("outPrefix (-o) must be specified")

des.table = sprintf("%s.txt", outPrefix)
des.plot = sprintf("%s.png", outPrefix)


tmp.df = NULL
nameL = NULL
for( src in srcL ){
	# sample=sampleL[1]
	name = sub(".spikeCnt.txt$", "", basename(src))
	nameL = c(nameL, name)
	tmp = read.delim(src, header=TRUE, row.names=1)

	if(is.null(tmp.df)){
		tmp.df = tmp[1:2,]
	}else{
		tmp.df = rbind(tmp.df, tmp[1:2,])
	}
}

#df = as.data.frame(t(tmp.df))
colnames(tmp.df) = c("Main","Spikein")
df = data.frame(Sample = nameL, tmp.df)


rownames(df) = NULL
df$Sample = factor(df$Sample, levels=nameL)
df$MainToSpikeRatio = df$Main / df$Spikein
df$SpikeFraction = df$Spikein / (df$Main + df$Spikein) * 100
df$ScaleFactor = mean(df$Spikein) / df$Spikein

write.table(df, des.table, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


g.Main = ggplot(data=df, aes(x=Sample, y=Main)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
	labs(title="Main", x="", y="Read depth")


g.Spike = ggplot(data=df, aes(x=Sample, y=Spikein)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
	labs(title="Spikein", x="", y="Read depth")

g.frac = ggplot(data=df, aes(x=Sample, y=SpikeFraction)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),) +
	labs(title="Spikein Fraction", x="", y="Spikein Fraction (%)")

g.ratio = ggplot(data=df, aes(x=Sample, y=MainToSpikeRatio)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),) +
	labs(title="Main-to-Spikein Ratio", x="", y="Reads Ratio: Main / Spike")

g.scale = ggplot(data=df, aes(x=Sample, y=ScaleFactor)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),) +
	labs(title="Scaling Factor (to multiply)", x="", y="Scaling Factor (to Multiply)")


g1 = plot_grid(g.Main, g.Spike, g.frac, ncol=1)
g2 = plot_grid(g.ratio, g.scale, ncol=1)
g = plot_grid(g1, g2, ncol=2)

ggsave(des.plot, g, width=3+0.3*nrow(df), height=8, dpi=200)




