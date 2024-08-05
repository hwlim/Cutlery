#!/usr/bin/env Rscript

#Create pdf/html report using R
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('plotly', quiet=TRUE))
suppressPackageStartupMessages(library('kableExtra', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))
suppressPackageStartupMessages(library('magick', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))

option_list <- list(
	make_option(c("-o","--outputFile"), help="prefix to output file; Can include path as well"),
	make_option(c("-g","--groupName"), help="Group name for the sample"),
	make_option(c("-s","--sampleDir"), help="Path to the Sample Folder, defined by 'sampleDir'"),
	make_option(c("-q","--qcDir"), help="Path to the Quality Control Folder, defined by 'qcDir'"),
	make_option(c("-f","--fragMixPrefix"), help="prefix for fragMix.txt file; output of frag_QC rule", default = "fragMix"),
	make_option(c("-t","--sampleTsvIn"), help="Path to sample.tsv")
)
parser <- OptionParser(usage = "%prog [options]",
	description="Description:
	Creates sample HTML report.
Input:
	'2.1.QualityControl' directory path
	Individual sample directory path
	Output file name
Output:
    - Report.html",
	option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)

# Option handling
opt=arguments$options
sampDir=normalizePath(opt$sampleDir)
qcDir=normalizePath(opt$qcDir)
groupName=opt$groupName
outputFile=opt$outputFile
sampleTsvIn=opt$sampleTsvIn
fragMixPrefix=opt$fragMixPrefix

#Read in and crop logo
logo = paste0(Sys.getenv("CUTLERY"), "/Plan/logo.png")
logo = image_read(logo)
logo = image_scale(logo, "220")
tempLogo = tempfile()
image_write(logo, path = tempLogo, format = "png")

#count peaks
countPeaks=function(bed) {
	countLines <- system(paste0("wc -l ", bed), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
	return(strtoi(strsplit(countLines, "\\s+")[[1]][1]))
}

#get sum of peak widths
sumPeaks=function(bed) {
	getDiffs <- system(paste0("awk '{$1 = $3 - $2} 1 {print $1}' ", bed), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
	return(sum(strtoi(getDiffs)))
}

#make two tables, one each for histone mode and factor mode
#make a fragQC table
histoneTable <- data.table(peakcount=numeric(), widths=numeric(), fraction=numeric())
factorTable <- data.table(peakcount=numeric(), fraction=numeric())
fragQCtable <- data.table(mode=character(), nfrFrac=numeric(), nucFrac=numeric(), w0=numeric(), w1=numeric(), w2=numeric())

#make a copy of the original Rmd file
template <- paste0(Sys.getenv("CUTLERY"), "/Script/cnr.createSampleReportHTMLTemplate.Rmd")
tempTemplate <- tempfile(fileext=".Rmd")
system(sprintf("cp %s %s", template, tempTemplate))

# Read sample.tsv
sampleIn = read.table(file = sampleTsvIn, sep = '\t', header = TRUE)

# Get sample name and sample QC directory
sampleName = tail(unlist(strsplit(sampDir, "/")), n=1)
sampleQC = paste0(sampDir, "/QC")

## get sample info from sample.tsv file
rowNum = which(sampleIn$Name == sampleName)
rowData = sampleIn[rowNum,]
peakMode = rowData$PeakMode

## get homer dir
homerFolderPath <- paste0(sampDir, "/", "HomerPeak.", peakMode)

#get fragmentQC files for the sample
fragQCfile <- paste0(sampleQC, "/", opt$fragMixPrefix ,".txt")
fragQC <- read.table(fragQCfile, header = TRUE)
fragLenDist <- paste0(sampleQC, "/fragLen.dist.png")

#get base frequency png file
baseFreq <- paste0(sampleQC, "/base_freq.png")

# skip visualization of control samples, but still retrieve fragQC
if (peakMode == "NULL") {
	peakMode = "ctrl"
	fragQCtable <- rbind(fragQCtable, list(peakMode, formatC(as.numeric(fragQC[[1]]), format="f", digits=2), 
	formatC(as.numeric(fragQC[[2]]), format="f", digits=2), formatC(as.numeric(fragQC[[3]]), format="f", digits=2),
	formatC(as.numeric(fragQC[[4]]), format="f", digits=2), formatC(as.numeric(fragQC[[5]]), format="f", digits=2)))
	next
}

if (peakMode == "factor") {
	heatmap <- paste0(homerFolderPath, "/heatmap.exBL.1rpm.png")
	peakFile <- paste0(homerFolderPath, "/peak.exBL.1rpm.bed")
	percFile <- paste0(homerFolderPath, "/peak.exBL.1rpm.stat")
	intersectPerc <- read.table(file = percFile, sep = "\t", header = FALSE)[[2]]
	factorTable <- rbind(factorTable, list(countPeaks(peakFile), paste0(intersectPerc, "%")))
} else {
	heatmap <- paste0(homerFolderPath, "/heatmap.exBL.png")
	peakFile <- paste0(homerFolderPath, "/peak.exBL.bed")
	percFile <- paste0(homerFolderPath, "/peak.exBL.stat")
	intersectPerc <- read.table(file = percFile, sep = "\t", header = FALSE)[[2]]
	histoneTable <- rbind(histoneTable, list(countPeaks(peakFile), sumPeaks(peakFile), paste0(intersectPerc, "%")))
}

#update fragQC table
fragQCtable <- rbind(fragQCtable, list(peakMode, formatC(as.numeric(fragQC[[1]]), format="f", digits=2), 
formatC(as.numeric(fragQC[[2]]), format="f", digits=2), formatC(as.numeric(fragQC[[3]]), format="f", digits=2),
formatC(as.numeric(fragQC[[4]]), format="f", digits=2), formatC(as.numeric(fragQC[[5]]), format="f", digits=2)))

#Render report
rmarkdown::render(tempTemplate, output_file = paste0(outputFile, ".html"), output_dir = sampleQC)

#Delete temp files
unlink(tempLogo)
unlink(tempTemplate)