#!/usr/bin/env Rscript

#Create pdf/html report using R
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('plotly', quiet=TRUE))
suppressPackageStartupMessages(library('kableExtra', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))
suppressPackageStartupMessages(library('magick', quiet=TRUE))

option_list <- list(
	make_option(c("-o","--outputFile"), help="prefix to output file; Can include path as well"),
	make_option(c("-s","--sampleDir"), help="Path to the Sample Folder, defined by 'sampleDir'"),
	make_option(c("-q","--qcDir"), help="Path to the Quality Control Folder, defined by 'qcDir'")
)
parser <- OptionParser(usage = "%prog [options]",
	description="Description:
	Creates final HTML report.
Input:
	'3.Sample' directory path
Output:
    - Report.html",
	 option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)

# Option handling
opt=arguments$options
sampDir=opt$sampleDir
qcDir=opt$qcDir
outputFile=opt$outputFile

#Read in and crop logo
logo = paste0(Sys.getenv("CUTLERY"), "/Plan/logo.png")
logo = image_read(logo)
logo = image_scale(logo, "200")
tempLogo = tempfile()
# logo = image_crop(logo, "240x100+33+10")
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
histoneTable <- data.table(Sample=character(), peakcount=numeric(), widths=numeric())
factorTable <- data.table(Sample=character(), peakcount=numeric())

#make a copy of the original Rmd file
template <- paste0(Sys.getenv("CUTLERY"), "/cnr.createReportHTMLTemplate.Rmd")
system(sprintf("cp %s Report.Rmd", template))

#create separate coverage html files for each sample
samplePaths <- list.files(sampDir)
sampleQC <- paste0(sampDir, "/", samplePaths, "/QC")

for (sample in sampleQC) {

	#get sample name
	sampleName <- tail(unlist(strsplit(sample, "/")), n=2)[1]

	#get path to sample Directory
	sampleDirectory <- paste(head(unlist(strsplit(sample, "/")), -1), collapse = "/")
	fileList <- list.files(sampleDirectory)
	
	#get heatmap and coverage
	homerFolderName <- grep("Homer", fileList, value=TRUE)
	factorMode <- tail(unlist(strsplit(homerFolderName, ".", fixed = TRUE)))[2]
	homerFolderPath <- paste0(sampleDirectory, "/", homerFolderName)
	
	#skip visualization of control samples
	if (is.null(factorMode)) {
		next
	}

	if (factorMode == "factor") {
		heatmap <- list.files(homerFolderPath, pattern = "*heatmap.exBL.1rpm.png", full.names = TRUE)
		peakFile <- list.files(homerFolderPath, pattern = "*peak.exBL.1rpm.bed", full.names = TRUE)
		factorTable <- rbind(factorTable, list(sampleName, countPeaks(peakFile)))
	} else {
		heatmap <- list.files(homerFolderPath, pattern = "*heatmap.exBL.png", full.names = TRUE)
		peakFile <- list.files(homerFolderPath, pattern = "*peak.exBL.bed", full.names = TRUE)
		histoneTable <- rbind(histoneTable, list(sampleName, countPeaks(peakFile), sumPeaks(peakFile)))
	}

	#get peak-examples png file
	peakExamplePlot <- list.files(homerFolderPath, pattern = "*peak.examples.png", full.names = TRUE)

	#write peak-examples section of Rmd file
	#add heatmap at the end
	#This is required as tabs need to be generated based on the number of samples in the sample.tsv file
	line = paste0("### ", sampleName, " {.tabset}")
	write(line,file="Report.Rmd",append=TRUE)
	line="```{r, echo=FALSE, out.width='100%', fig.align='center', message=FALSE, warning=FALSE}"
	write(line,file="Report.Rmd",append=TRUE)
	line = paste0("knitr::include_graphics(\"",peakExamplePlot,"\")")
	write(line,file="Report.Rmd",append=TRUE)
	line = "```"
	write(line,file="Report.Rmd",append=TRUE)

	line = "<br />"
	write(line,file="Report.Rmd",append=TRUE)

	line="```{r, echo=FALSE, out.width='70%', fig.align=\"center\"}"
	write(line,file="Report.Rmd",append=TRUE)
	line = paste0("knitr::include_graphics(\"",heatmap,"\")")
	write(line,file="Report.Rmd",append=TRUE)
	line = "```"
	write(line,file="Report.Rmd",append=TRUE)

	line = "\n"
	write(line,file="Report.Rmd",append=TRUE)
	line = "\n"
	write(line,file="Report.Rmd",append=TRUE)
	line = "<br />"
	write(line,file="Report.Rmd",append=TRUE)
	line = "\n"
	write(line,file="Report.Rmd",append=TRUE)
	line = "\n"
	write(line,file="Report.Rmd",append=TRUE)
}

#create master Report Template
#includes everything except for genome coverage plots
rmarkdown::render('Report.Rmd', output_file = paste0(outputFile, ".html"))

#Delete intermediate coverage Reports and plots
system(sprintf("rm Report.Rmd"))

unlink(tempLogo)
