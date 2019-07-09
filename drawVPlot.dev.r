#!/usr/bin/env Rscript


############################################
# Written by Hee-Woong Lim
#
# Draw V-plot(s) using given anchor bed file and center/fragLen file
#
# Things to consider:
# - When to select fragments overlapping with anchor regions
# 	In this script? or in creating cfl files? or both?

suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('dplyr', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))

source("~/bin/commonR.r")
source("~/bin/basicR.r")
source("~/bin/scRNAseq/seurat.common.r")

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="seurat.reduce", help="output prefix potentially including path"),
#	make_option(c("-M","--maxDim"), default=30, help="Max number of dimensions for PCA calculation, must be >= N.dim. [default: 30]"),
	make_option(c("-N","--N.dim"), default=20, help="Number of dimensions to use from PCA for tSNE/UMAP, must be <= maxDim. [default: 20]"),
	make_option(c("-n","--N.embed"), default=2, help="Number of dimensions for tSNE/UMAP, must be 2 or 3. [defalt: 2]"),
	make_option(c("-P","--pcaOnly"), default=FALSE, action="store_true", help="Run PCA only. In default, All of PCA/tSNE/UMAP are performed. [default: Off]"),
	make_option(c("-s","--skipPCA"), default=FALSE, action="store_true", help="Skip PCA and do tSNE/UMAP only. [default: Off]"),
	make_option(c("-d","--downSampleBy"), default=NULL, help="Down sampling criteria. NULL / orig.ident / group. If set, Down sampled version of DimPlot will be created. [default: NULL]"),
	make_option(c("-c","--cells"), default=NULL, help="Plain text file containing cell list to select before dimension reduction. [default: NULL]")
#	make_option(c("-t","--thread"), default=1, help="Number of thread for JackStraw [default: 1]")
#	make_option(c("-o","--outfile"), default="barplot.png", help="Output file, with .png extension. [default: barplot.png]")
#	make_option(c("-f","--field"), default="", help="Comma-separated field numbers for x-axis, y-axis."),
)
parser <- OptionParser(usage = "%prog [options] <seurat R object file>",
			description="Performs cell clustering in Seurat and create following results:\
	- <outPrefix>.rds : Seurat object file. possibly overwrite the source rds file\
	- <outPrefix>.1.1.PCA.png
	- <outPrefix>.1.2.JackStrawPlot.png",
			option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires a data file")
} else {
	src <- arguments$args[1]
}

###########################
## Option handling
opt=arguments$options
outPrefix=opt$outPrefix
maxDim = opt$maxDim
N.dim = opt$N.dim
dim.embed = opt$N.embed
pcaOnly = opt$pcaOnly
skipPCA = opt$skipPCA
#thread = opt$thread
downSampleBy = opt$downSampleBy
src.cells = opt$cells

if(FALSE){
	outPrefix="test/seurat"
	maxDim=30
	N.dim = 10
	dim.embed=2
	pcaOnly=FALSE
	skipPCA=FALSE
#	thread = 4
	downSampleBy = NULL
	src.cells=NULL
	src = "seurat.rds"
}
###########################
## Option validation
stopifnot(maxDim >= dim.embed)
stopifnot(maxDim > 0 && dim.embed > 0)
stopifnot( is.null(downSampleBy) || downSampleBy %in% c("orig.ident", "group") )
stopifnot( ! (pcaOnly && skipPCA) )
assertFileExist(src)
if(!is.null(src.cells)) assertFileExist(src.cells)
des.seurat = sprintf("%s.rds", outPrefix)
desDir=dirname(outPrefix)





write(sprintf("=============================="), stderr())
write(sprintf("Performing dimension reduction"), stderr())
write(sprintf("  maxDim = %s", maxDim), stderr())
write(sprintf("  N.dim = %s", N.dim), stderr())
write(sprintf("  dim.embed = %s", dim.embed), stderr())
write(sprintf("  pcaOnly = %s", pcaOnly), stderr())
write(sprintf("  skipPCA = %s", skipPCA), stderr())