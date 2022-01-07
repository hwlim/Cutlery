#!/usr/bin/env Rscript

# source("/Volumes/limlab/ChristopherAhn/21.09.24_create_pdf_py/genomeR.r")
# source("/Volumes/limlab/ChristopherAhn/21.09.24_create_pdf_py/basicR.r")
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('reshape2', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))

option_list <- list(
    make_option(c("-o","--outPrefix"), help="Path + prefix for output file"),
    make_option(c("-m","--mode"), help="Peak mode; Enter 'histone' or 'factor'. This will modify the bin size and window size used when extracting data from bw files"),
    make_option(c("-n","--numPeaks"), help="Number of highest scoring peaks to visualize; should be an integer."),
    make_option(c("-p","--peak"), help="Path to peak.bed file")
)

parser <- OptionParser(usage = "%prog",
	description="Description:
	Takes a peak file to visualize the top scoring peak regions and save as one png and one pdf file.
Input:
	sample/directory
Output:
    - <snapshot>.png",
	 option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)

# Option handling
opt=arguments$options
outPrefix=opt$outPrefix
mode=opt$mode
numHighestPeaks=opt$numPeaks
peakFile <- opt$peak
bwFiles <- arguments$args


#Create NFR and NUC plots and merge
create_plots_and_add_to_plot_list=function(tempNUC.melt, tempNFR.melt, nucMax, nfrMax, col, chr){
    
    #check if max RPM is zero in both the sample and the control. If so, ylim has to be manually set
    nucMax <- max(tempNUC.melt$value, na.rm = TRUE)
    nfrMax <- max(tempNFR.melt$value, na.rm = TRUE)

    #if max RPM in nfr bw is zero, set ylim manually
    NFR <- draw_NFR(tempNFR.melt, nfrMax, col)
    
    #if max RPM in nuc bw is zero, set ylim manually
    NUC <- draw_NUC(tempNUC.melt, nucMax, chr)

    #merge both nfr and nuc visualizations vertically
    merged <- plot_grid(NFR, NUC, ncol = 1, align = "v")

    #create y-axis label
    newMerged <- ggdraw(merged) +
    draw_label("RPM", x=0.01, y=0.5, vjust=0.5, hjust=0.5, angle = 90, fontface="bold", size = 12)
    
    return(list(newMerged))
}

#Draw one NFR plot
draw_NFR=function(tempNFR.melt, nfrMax, col){
    
    NFR <- ggplot(tempNFR.melt, aes(rows,value)) +
        geom_area(colour=tempNFR.melt$color, fill=tempNFR.melt$color) + 
        facet_grid(series ~ .) + 
        labs(title=paste0("Highest Scoring Peak #", col), y = "Score") +
        theme(plot.title = element_text(hjust = 0.5, face="bold"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(colour="white")
            ) +
        theme(strip.background = element_rect(fill="red")) +
        theme(strip.text = element_text(colour = 'white', face="bold", size = 14)) +
        {if(nfrMax == 0)ylim(0,5)}
    
    return(NFR)
}

#Draw one NUC plot
draw_NUC=function(tempNUC.melt, nucMax, chr){
    
    NUC <- ggplot(tempNUC.melt, aes(rows,value)) +
        geom_area(colour=tempNUC.melt$color, fill=tempNUC.melt$color) + 
        facet_grid(series ~ .) + 
        labs(x = chr, y = "Score") +
        theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill=NA),
            axis.title.x = element_text(face="bold"),
            axis.title.y = element_text(colour="white")
            ) +
        theme(strip.background = element_rect(fill="blue")) +
        theme(strip.text = element_text(colour = 'white', face="bold", size = 14)) +
        {if(nucMax == 0)ylim(0,5)}
    
    return(NUC)
}

draw_peak_example=function(bed, extractedBW, windowSize, binsize, des=NULL){
#draw coverage plot for all peaks in a bed file. Requires a bed file and at least two bigwig files.
#ctrl bigwig files are "NULL" by default, unless specified
#binsize for extractBigWigData is 20 by default
#window size for extractBigWigData is 10000 by default
#required packages: ggplot2, reshape2, cowplot

    #If there is a ctrl for this sample
    if (length(extractedBW) == 4){
        #convert to DF and transpose
        nfrdf <- as.data.frame(t(extractedBW[[1]]))
        nucdf <- as.data.frame(t(extractedBW[[2]]))
        ctrlNfrdf <- as.data.frame(t(extractedBW[[3]]))
        ctrlNucdf <- as.data.frame(t(extractedBW[[4]]))
        tempNFR <- as.data.frame(t(extractedBW[[1]]))
        tempNUC <- as.data.frame(t(extractedBW[[2]]))

        #ratio is used for final plot aspect ratio
        ratio = 6

        #Add x axis numbering
        tempNUC <- cbind(rows = 1:nrow(tempNUC), tempNUC)

        #extract first three cols
        #this is due to errors when trying to create a new dataframe. Tedious work...
        tempNUC = tempNUC[,1:3]

        #rename column headers
        colnames(tempNUC) <- c("rows", "NUC", "NUC_Ctrl")

        #Add x axis numbering
        tempNFR <- cbind(rows = 1:nrow(tempNFR), tempNFR)
        
        #extract first three cols
        #this is due to errors when trying to create a new dataframe. Tedious work...
        tempNFR = tempNFR[,1:3]

        #rename column headers
        colnames(tempNFR) <- c("rows", "NFR", "NFR_Ctrl")

    #if there is no control for this particular sample
    } else {
        nfrdf <- as.data.frame(t(extractedBW[[1]]))
        nucdf <- as.data.frame(t(extractedBW[[2]]))
        tempNFR <- as.data.frame(t(extractedBW[[1]]))
        tempNUC <- as.data.frame(t(extractedBW[[2]]))

        #ratio is used for final plot aspect ratio
        ratio = 3

        #Add x axis numbering
		tempNUC <- cbind(rows = 1:nrow(tempNUC), tempNUC)
		
        #extract first two cols
		#this is due to errors when trying to create a new dataframe. Tedious work...
        tempNUC = tempNUC[,1:2]

		#rename column headers
		colnames(tempNUC) <- c("rows", "NUC")

		#Add x axis numbering
		tempNFR <- cbind(rows = 1:nrow(tempNFR), tempNFR)
		
        #extract first two cols
		#this is due to errors when trying to create a new dataframe. Tedious work...
        tempNFR = tempNFR[,1:2]

		#rename column headers
		colnames(tempNFR) <- c("rows", "NFR")
    }

    #read peak file
    peakFile <- read.table(bed, header = FALSE, sep = "\t")

    #create empty list of images
    listOfPlots <- list()

    #append images vertically
    #create plots for each highest scoring peak
    for (col in 1:ncol(nucdf)) {

        #get peak start
        peakStart <- peakFile[col,2]
        chr <- peakFile[col,1]

        #get x axis range
        xRange <- peakStart:(peakStart + (windowSize/binsize) - 1)

        #modify column values
        tempNUC[1] <- xRange * binsize
        tempNFR[1] <- xRange * binsize

        #if control exists
        if (length(extractedBW) == 4){
            tempNFR[2] <- nfrdf[col]
            tempNFR[3] <- ctrlNfrdf[col]
            tempNUC[2] <- nucdf[col]
            tempNUC[3] <- ctrlNucdf[col]

        #if control doesn't exist
        } else {
            tempNFR[2] <- nfrdf[col]
			tempNUC[2] <- nucdf[col]
        }

        #replace zeros to NA
        tempNFR[tempNFR == NA] <- 0
        tempNUC[tempNUC == NA] <- 0

        tempNUC.melt <- melt(tempNUC, id.vars = 'rows', variable.name = 'series')
        tempNFR.melt <- melt(tempNFR, id.vars = 'rows', variable.name = 'series')

        #if control exists
        if (length(extractedBW) == 4){
            #define an array of type for color coding; this is for samples with controls
            colorTypeNUC=rep(c("blue", "grey43"), each=(nrow(tempNUC.melt)/2))
            colorTypeNFR=rep(c("red", "grey43"), each=(nrow(tempNFR.melt)/2))

        } else {
            #define an array of type for color coding; this is for samples w/o controls
            colorTypeNUC=rep(c("blue"), each=(nrow(tempNUC.melt)))
			colorTypeNFR=rep(c("red"), each=(nrow(tempNFR.melt)))
        }

        #add colors to df
        tempNUC.melt$color <- colorTypeNUC
        tempNFR.melt$color <- colorTypeNFR

        #create plots for one peak by drawing and merging the NUC and NFR plots
        #append current sample's peak example visualization to a list
        listOfPlots <- append(listOfPlots, create_plots_and_add_to_plot_list(tempNUC.melt, tempNFR.melt, nucMax, nfrMax, col, chr))

    } #end of loop to create plots

    #merge all peak examples
    allMerged <- plot_grid(plotlist=listOfPlots, ncol = 1, align = "v")

    #save peak examples as png and pdf
    ggsave(paste0(des, ".png"), allMerged, width = 10, height = (strtoi(numHighestPeaks) * ratio), dpi = 500, units = "in", device = "png")
    ggsave(paste0(des, ".pdf"), allMerged, width = 10, height = (strtoi(numHighestPeaks) * ratio), dpi = 500, units = "in", device = "pdf")

}

#__main__
sampTemp <- head(unlist(strsplit(peakFile, "/")), -2)
sampDir <- paste(sampTemp, collapse = '/')

#Parse string to get sample Name
sampleName <- tail(unlist(strsplit(sampDir, "/")), n=1)

#Set temp peak.bed file for top peaks
tmpBed = tempfile()
system(sprintf(paste0("head -n ", numHighestPeaks, " ", peakFile, " > ", tmpBed)))

#Set window size and bin size for bw extraction
if (mode == 'factor') {
    windowSize = 1000
    binsize = 2
} else {
    windowSize = 10000
    binsize = 20
}

#extractBW
extractedBW = list()
for (i in seq_along(bwFiles)) {
    extractedBW[[i]] = extractBigWigData(bed = tmpBed, bw = bwFiles[i], width = windowSize, binSize = binsize)
}

#Draw peak examples
draw_peak_example(tmpBed, extractedBW, windowSize = windowSize, binsize = binsize, des=outPrefix)

#remove temp peak.bed file
unlink(tmpBed)