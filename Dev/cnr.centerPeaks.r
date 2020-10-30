#!/usr/bin/env Rscript

#### NOTE
## After extracting profile and visual inspection,
## I noticed that most of the peak are well-centered not only at the maximum stack-up but also in the middle of maximum plateau
## except for the situation of more than one plateau with the same maximum height
##




suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressMessages(suppressPackageStartupMessages(library('Biostrings', quiet=TRUE)))
suppressMessages(suppressPackageStartupMessages(library('tools', quiet=TRUE)))


source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))

option_list <- list(
	make_option(c("-o","--out"), default=NULL, help="Output files, default=<input file>.ctr.bed"),
	make_option(c("-n","--n_wma"), default=0, help="Bandwidth for WMA, half-window size, default=0")
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
	- Centered peak file containing the same sized regions" )

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
des=opt$out
n.wma=opt$n_wma

stopifnot(is.numeric(n.wma) && n.wma>=0)

if(FALSE){
	src.peak="../TestData/HNF4a_Park/peak.exBL.1rpm.bed"
	#src.bw1="../TestData/HNF4a_Park/igv.nfr.ctr.bw"
	#src.bw2="../TestData/HNF4a_Park/igv.allFrag.bw"
	src.bw="../TestData/HNF4a_Park/igv.nfr.con.bw"
	des="test.centered.bed"
	n.wma=5
}

assertFileExist(c(src.peak, src.bw))
if( is.null(des) ) des = sprintf("%s.ctr", sub(".bed$","", basename(src.peak)))


write(sprintf("Centering CUT&RUN peaks"), stderr())
write(sprintf("  - src = %s", src.peak), stderr())
write(sprintf("  - bigwig = %s", src.bw), stderr())
write(sprintf("  - bandwidth for WMA = %d", n.wma), stderr())
write(sprintf("  - des = %s", des), stderr())


write(sprintf("Reading data..."), stderr())
peak = readBedFile(src.peak)
## Homer peak is sometimes incorrect sized
## Forced resizing to handle those incorrect sized peaks
peak = resizeBed(peak,200)
half_width = round(peak[1,3] - peak[1,2])/2

data = extractBigWigData1bp(peak, src.bw, des=NULL)
stopifnot(sum(is.na(data))==0)

#stopifnot(all(peak[,3]-peak[,2] == peak[1,3]-peak[1,2]))
if(FALSE){
	## Dev. code to compare different types of bigwig files
	data1 = extractBigWigData1bp(src.peak, src.bw1, des=NULL)
	data2 = extractBigWigData1bp(src.peak, src.bw2, des=NULL)
	data3 = extractBigWigData1bp(src.peak, src.bw3, des=NULL)

	draw=function(i){
		y1 = as.numeric(data1[i,])
		y1 = y1/max(y1)
		y2 = as.numeric(data2[i,])
		y2 = y2/max(y2)
		y3 = as.numeric(data3[i,])
		y3 = y3/max(y3)

		tmp = data.frame(nfr_ctr = y1, all = y2, nfr_con=y3 )
		matplot(tmp, type="l")
	}
	draw(1)

	testWMA=function(i=0){
		if(i==0){
			y = as.numeric(data[sample(nrow(data),1),])
		}else {
			y = as.numeric(data[i,])
		}
		tmp = data.frame(y, wma(y,5), wma(y,10), wma(y,20))
		matplot(tmp,type="l")
	}
}


findCenter=function(x, n.wma=0){
	## n.wma
	if(n.wma==0){
		ctr = mean(which(x==max(x)))
	}else{
		x = wma(x, n.wma)
		ctr = mean(which(x==max(x)))
	}
}

## NOTE: Some center after WMA goes beyond the initial peak boundary
write(sprintf("Centering..."), stderr())
offset = apply(data,1, function(x) findCenter(x,n.wma))
offset = round(offset)
#hist(center,100)
#tmp.index=which(center==200)

center = peak[,2] + offset
peak[,2] = center - half_width
peak[,3] = center + half_width


write(sprintf("Writing..."), stderr())
writeBedFile(peak, des)
