#!/usr/bin/env Rscript

#suppressPackageStartupMessages(library('geneplotter', quiet=TRUE))
#suppressPackageStartupMessages(library('RColorBrewer', quiet=TRUE))
#suppressPackageStartupMessages(library('MASS', quiet=TRUE))
#suppressPackageStartupMessages(library('robustbase', quiet=TRUE))
#suppressPackageStartupMessages(library('KernSmooth', quiet=TRUE))
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
#suppressMessages(suppressPackageStartupMessages(library('Biostrings', quiet=TRUE)))
suppressMessages(suppressPackageStartupMessages(library('tools', quiet=TRUE)))


source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/motifR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/IDOM/IDOM.r", Sys.getenv("LIMLAB_BASE")))
source(sprintf("%s/ExoTools/ChipExoUtil.r", Sys.getenv("LIMLAB_BASE")))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default=NULL, help="Prefix of output files, default=< anchor file name except for path and extension"),
	make_option(c("-m","--margin"), default=40, help="Flanking window size to consider (bp). default=40"),
	make_option(c("-p","--pseudo"), default=0.1, help="Pseudo value to add to calculate contrast, default=0.1"),
	make_option(c("-f","--offset"), default="0,0", help="Comma-separated offset to reduce/increase anchor regions. <left>,<right> (primarily to trim extra portion of motif), default=0,0")
#	make_option(c("-k","--keepData"), default=FALSE, action="store_true", help="Keep extracted data, default=FALSE"),
#	make_option(c("-v","--verbose"), default=FALSE, help="Verbose"),
)
parser <- OptionParser(usage = "%prog [options] [anchor bed file] [bigWig prefix]", option_list=option_list,
	description="Measure & compare footprint signal between anchor window vs flanking regions
Input:
	- anchor bed file
	- bigWig file prefix for a pair strands, i.e.<prefix>.<plus/minus>.bw
Output:
	- Text file of anchor bed concatenated with footprint signal at anchor adn flanks + contrast (anchor/flank)
	[anchor] / average within anchor & flanks / Contrast=log2( (flank + pseudo)/(anchor + pseudo) )
	- Unified anchors (unique 4th column names) by selecting maximum contrast" )

#################################
# Option handling
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) < 2){
	write(sprintf("Error: Requires at bed and bigWig prefix"), stderr())
	print_help(parser)
	q()
}else{
	src.anchor=arguments$args[1]
	src.bwPrefix = arguments$args[2]
}
opt=arguments$opt
outPrefix=opt$outPrefix
margin=opt$margin
pseudo=opt$pseudo
offset=opt$offset

if(FALSE){
	outPrefix=NULL
	src.anchor="../GSX2_D_1.300bp.p5e-4.bed"
	src.bwPrefix="../../../E12.5_Flag.1bp"
	pseudo=0.1
	offset="1,3"
	margin=40

}


offset=as.numeric(strsplit(offset,",")[[1]])
stopifnot(all(is.numeric(offset)))

src.bwPlus = sprintf("%s.plus.bw", src.bwPrefix)
src.bwMinus = sprintf("%s.minus.bw", src.bwPrefix)
assertFileExist(c(src.anchor, src.bwPlus, src.bwMinus))



if( is.null(outPrefix) ) outPrefix = sub(".bed$","", basename(src.anchor))
outDir=dirname(outPrefix)
system(sprintf("mkdir -p %s", outDir))



## DNA sequence
write(sprintf("================================="), stderr())
write(sprintf("Footprint contrast"), stderr())
write(sprintf("================================="), stderr())
write(sprintf("anchor = %s", src.anchor), stderr())
write(sprintf("bwPrefix = %s", src.bwPrefix), stderr())
write(sprintf("margin = %d", margin), stderr())
write(sprintf("pseudo = %s", pseudo), stderr())
write(sprintf("offset = %s", paste(offset, collapse=",")), stderr())
write(sprintf("outPrefix = %s", outPrefix), stderr())


des.bwPlus =  sprintf("%s.profile.plus.gz", outPrefix)
des.bwMinus = sprintf("%s.profile.minus.gz", outPrefix)


data.anchor=readBedFile(src.anchor)[,1:6]
colnames(data.anchor) = c("Chr","Start","End","Name","Null","Direc")

anchor.resized = data.anchor
tf.minus = anchor.resized[,6]=="-"
anchor.resized[!tf.minus,2] = anchor.resized[!tf.minus,2] + offset[1]
anchor.resized[!tf.minus,3] = anchor.resized[!tf.minus,3] - offset[2]
anchor.resized[tf.minus,2] = anchor.resized[tf.minus,2] + offset[2]
anchor.resized[tf.minus,3] = anchor.resized[tf.minus,3] - offset[1]

anchorSize = anchor.resized[1,3]-anchor.resized[1,2]
stopifnot(all(anchorSize == anchor.resized[,3]-anchor.resized[,2]))
anchor.uniqName = data.frame( anchor.resized[,1:3], sprintf("anchor.%d", 1:nrow(anchor.resized)), anchor.resized[,c(5,6,4)] )



#####################################################3
## BigWig footprint
## Note:
##	Initially, bigWig signal is extracted in unstranded manner.
##	They are switched and flipped manually according to the direction in the next step. 
write(sprintf("2) Extracting footprint signal"), stderr())
#if( !all(file.exists( des.bwPlus, des.bwMinus )) ){
bed.profile = tempfile()
bed.ext = anchor.uniqName
bed.ext[,2] = bed.ext[,2] - margin
bed.ext[,3] = bed.ext[,3] + margin
bed.ext[,5] = 0
bed.ext[,6] = "+"
write.table(bed.ext[,1:6], bed.profile, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
	
write(sprintf("\tplus = %s", des.bwPlus), stderr())
cmd = sprintf("bwtool extract bed -fill=0 -decimals=5 -tabs %s %s /dev/stdout | cut -f 4,8- | gzip -c > %s", bed.profile, src.bwPlus, des.bwPlus)
system(cmd)
write(sprintf("\tminus = %s", des.bwMinus), stderr())
cmd = sprintf("bwtool extract bed -fill=0 -decimals=5 -tabs %s %s /dev/stdout | cut -f 4,8- | gzip -c > %s", bed.profile, src.bwMinus, des.bwMinus)
system(cmd)
unlink(bed.profile)

data.profileL = readExoProfile2(srcPlus=des.bwPlus, srcMinus=des.bwMinus, direc=data.anchor[,6])

idx.anchor = (margin + 1):(margin+anchorSize)
idx.flank = c( 1:margin, (margin+anchorSize+1):(2*margin+anchorSize) )

signal.anchor = apply( data.profileL[["plus"]][,idx.anchor], 1, mean ) + abs(apply( data.profileL[["minus"]][,idx.anchor], 1, mean ))
signal.flank = apply( data.profileL[["plus"]][,idx.flank], 1, mean ) + abs(apply( data.profileL[["minus"]][,idx.flank], 1, mean ))


### All data
data.final = data.frame(data.anchor,
				Anchor=signal.anchor,
				Flank = signal.flank,
				Contrast = log2( (signal.flank + pseudo) / (signal.anchor + pseudo) )
			)
write.table(data.final, sprintf("%s.all.txt", outPrefix), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

### 2D density plot
tmp = log2(data.final[,c("Flank","Anchor")] + pseudo)
axisMin=min(tmp)
axisMax=max(tmp)
axisLim=c(axisMin, axisMax)

png(sprintf("%s.density.png", outPrefix), width=500, height=550)
par(las=1, cex=1.5)
drawDensity2D(log2(data.final[,c("Flank","Anchor")] + pseudo), xlim=axisLim, ylim=axisLim)
title(xlab=sprintf("Flanking average signal (%d bp)", margin),
	ylab=sprintf("Anchor average signal (%d bp)", anchorSize),
	main="Footprint Contrast")
dev.off()
