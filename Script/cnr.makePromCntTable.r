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
	make_option(c("-o","--outPrefix"), default=NULL, help="Output prefix. required")
#	make_option(c("-f","--bamFlag"), default="0x2", help="flag for bam records. NULL is allowed to unset. Ignored for bed file. default=0x2 (concordant pairs only)"),
#	make_option(c("-F","--bamUnFlag"), default="0x400", help="flag for bam records to exclude, NULL is allowed to unset. Ignored for bed file. default=0x400 (exclude duplicates)")
#	make_option(c("-t","--title"), default="Title", help="Main Title [default: Title]"),
#	make_option(c("-s","--size"), default="600,600", help="Comma-separated figure size, xSize,ySize"),
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <txt1> <txt2> ...",
	description="Description:
	Collect target tag count and spikein count information from multiple files and make a table file and visualize bar plots
	Assuming each input file is in the format:
		Name  ****
		TTC       ####
		PromoterCnt    ####
		PromoterFrac    ####
Output:
	<outPrefix>.txt: A tabular txt file summarizing total / unique fragment counts
	<outPrefix>.pdf: A series of diagnostic plots",
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
		"3.Sample/E12.5_Flag_1/QC/fragment.uniq_cnt.txt",
"3.Sample/E12.5_Flag_2/QC/fragment.uniq_cnt.txt",
"3.Sample/E12.5_H3K27me3_Flag_1/QC/fragment.uniq_cnt.txt",
"3.Sample/E12.5_H3K27me3_Flag_2/QC/fragment.uniq_cnt.txt",
"3.Sample/E12.5_H3K27me3_early/QC/fragment.uniq_cnt.txt",
"3.Sample/E12.5_IgG_1/QC/fragment.uniq_cnt.txt",
"3.Sample/E12.5_IgG_2/QC/fragment.uniq_cnt.txt",
"3.Sample/E13.5_Gsx2/QC/fragment.uniq_cnt.txt",
"3.Sample/E13.5_Gsx2_1to1000/QC/fragment.uniq_cnt.txt",
"3.Sample/E13.5_IgG/QC/fragment.uniq_cnt.txt"
	)
	outPrefix="test"
}
assertFileExist(srcL)
if(is.null(outPrefix)) stop("outPrefix (-o) must be specified")

desDir=dirname(outPrefix)
system(sprintf("mkdir -p %s", desDir))
des.txt = sprintf("%s.txt", outPrefix)
des.pdf = sprintf("%s.pdf", outPrefix)
des.png = sprintf("%s.png", outPrefix)

data = NULL
## 3 x N data.frame: row -> statistics; column -> sample
for( src in srcL ){
	# src=srcL[1]
	## 1 column data with sample name in the column name
	tmp = read.delim(src, header=TRUE, row.names=1, check.names=F)
	if(is.null(data)){
		data = tmp
	}else{
		data = cbind(data, tmp)
	}
}
data = t(data)
data = data.frame(Name=rownames(data), data)
data$Name = factor(data$Name, levels=data$Name)

write.table(data, des.txt, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


g.ttc = ggplot(data=data, aes(x=Name, y=TTC)) +
	geom_bar(stat="identity", fill="steelblue")+
	#theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
	labs(title="Total Fragments", x="", y="Count") +
	coord_flip()

g.uniq = ggplot(data=data, aes(x=Name, y=Promoter)) +
	geom_bar(stat="identity", fill="steelblue")+
	#theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
	labs(title="Unique Fragments", x="", y="Count") +
	coord_flip()

g.frac = ggplot(data=data, aes(x=Name, y=Percentage)) +
	geom_bar(stat="identity", fill="steelblue")+
	#theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),) +
	labs(title="% of Unique Fragments", x="", y="Unique Fragments (%)") +
	coord_flip()

g.scatter = ggplot(data=data, aes(x=log10(TTC), y=Percentage)) +
	geom_point() + ylim(0,100) + xlim(0,max(8, max(log10(data$TTC)))) +
	#theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),) +
	labs(title="TTC vs Promoter Fraction", x="log10(Total Fragments)", y="Promoter (%)")

g.scatter.labeled = g.scatter +
	geom_text(label=rownames(data), size=2) +
	labs(title="TTC vs UniqFrac, Labeled")


g1 = plot_grid(g.ttc, g.uniq, g.frac, ncol=1)
g2 = plot_grid(g.scatter, g.scatter.labeled, ncol=1)
g = plot_grid(g1, g2, ncol=2)


ggsave(des.pdf, g, width=12, height=12)
ggsave(des.png, g, width=12, height=12)



