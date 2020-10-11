#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressMessages(suppressPackageStartupMessages(library('Biostrings', quiet=TRUE)))
suppressMessages(suppressPackageStartupMessages(library('tools', quiet=TRUE)))


source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))

option_list <- list(
	make_option(c("-o","--outPrefix"), default=NULL, help="Prefix of output files, default=<input file>.ctr.bed")
#	make_option(c("-v","--verbose"), default=FALSE, help="Verbose"),
)
parser <- OptionParser(usage = "%prog [options] [peak] [bigwig]", option_list=option_list,
	description="Description:
	Find center of peaks with maximum fragment stack-up from a given bigwig file and create a new peak file
	If more than one max position, their middle (average) location will be the center
Input:
	- Peak file in bed format containing fixed size regions
	- BigWig file of fragment stack-height profile (i.e. igv.nfr.ctr.bw)
Output:
	- <outPrefix>.bed: Centered peak file containing the same sized regions" )

#################################
# Option handling
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) < 2){
	write(sprintf("Error: Requires at bed and bigWig prefix: (%d)", length(arguments$args)), stderr())
	print_help(parser)
	q()
}else{
	src.peak=arguments$args[1]
	src.bw = arguments$args[2]
}
opt=arguments$opt
outPrefix=opt$outPrefix


if(FALSE){
	src.peak="../TestData/HNF4a_Park/peak.exBL.1rpm.bed"
	src.bw="../TestData/HNF4a_Park/igv.nfr.ctr.bw"
	outPrefix=NULL
}

assertFileExist(c(src.peak, src.bw))
if( is.null(outPrefix) ) outPrefix = sprintf("%s.ctr", sub(".bed$","", basename(src.peak)))



data = extractBigWigData1bp(src.peak, src.bw, des=NULL)
## After extracting profile and visual inspection,
## I noticed that most of the peak are well-centered not only at the maximum stack-up but also in the middle of maximum plateau
## except for the situation of more than one plateau with the same maximum height
