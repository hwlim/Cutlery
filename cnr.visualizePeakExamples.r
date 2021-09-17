#!/usr/bin/env Rscript

# source("/Volumes/limlab/ChristopherAhn/create_pdf_py/genomeR.r")
# source("/Volumes/limlab/ChristopherAhn/create_pdf_py/basicR.r")
source(sprintf("%s/genomeR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('reshape2', quiet=TRUE))
suppressPackageStartupMessages(library('cowplot', quiet=TRUE))

option_list <- list(
    make_option(c("-o","--outDir"), help="Path to output directory"),
    make_option(c("-m","--mode"), help="Peak mode; Enter 'histone' or 'factor'"),
    make_option(c("-c","--ctrlName"), help="Name of control sample; 'NULL' if there is no control sample"),
    make_option(c("-n","--numPeaks"), help="Number of highest scoring peaks to visualize; should be an integer")
)

parser <- OptionParser(usage = "%prog",
	description="Description:
	Takes a peak file to visualize the top 6 high scoring peak regions and save as images.
Input:
	sample/directory
Output:
    - <snapshot>.png",
	 option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)

# Option handling
opt=arguments$options
outDir=opt$outDir
mode=opt$mode
ctrlSampleName=opt$ctrlName
numHighestPeaks=opt$numPeaks
peakFilePath <- arguments$args

draw_peak_example=function(bed, nucBW, nfrBW, ctrlNucBW, ctrlNfrBW, windowSize = 10000, binsize = 20, des=NULL){
#draw coverage plot for all peaks in a bed file. Requires a bed file and at least two bigwig files.
#ctrl bigwig files are "NULL" by default, unless specified
#binsize for extractBigWigData is 20 by default
#window size for extractBigWigData is 10000 by default
#required packages: ggplot2, reshape2, cowplot
    
    #if control sample exists
    if (!is.null(ctrlNucBW)) {

		#extract by 10bp
		nucPeakTable <- extractBigWigData(bed = bed, bw = nucBW, width = windowSize, binSize = binsize)
		nfrPeakTable <- extractBigWigData(bed = bed, bw = nfrBW, width = windowSize, binSize = binsize)
		nucCtrlTable <- extractBigWigData(bed = bed, bw = ctrlNucBW, width = windowSize, binSize = binsize)
		nfrCtrlTable <- extractBigWigData(bed = bed, bw = ctrlNfrBW, width = windowSize, binSize = binsize)

		#convert to DF and transpose
		nucdf <- as.data.frame(t(nucPeakTable))
		nfrdf <- as.data.frame(t(nfrPeakTable))
		ctrlNucdf <- as.data.frame(t(nucCtrlTable))
		ctrlNfrdf <- as.data.frame(t(nfrCtrlTable))
		tempNUC <- as.data.frame(t(nucPeakTable))
		tempNFR <- as.data.frame(t(nfrPeakTable))

		#Add x axis numbering
		tempNUC <- cbind(rows = 1:nrow(tempNUC), tempNUC)

		#remove unnecessary columns
		#this is due to errors when trying to create a new dataframe. Tedious work...
        tempNUC = tempNUC[,1:3]

		#rename column headers
		colnames(tempNUC) <- c("rows", "NUC", "NUC_Ctrl")

		#Add x axis numbering
		tempNFR <- cbind(rows = 1:nrow(tempNFR), tempNFR)
		#this is due to errors when trying to create a new dataframe. Tedious work...
		tempNFR = tempNFR[,1:3]

		#rename column headers
		colnames(tempNFR) <- c("rows", "NFR", "NFR_Ctrl")

		#read peak file
		peakFile <- read.table(bed, header = FALSE, sep = "\t")

        #create empty list of images
        listOfPlots <- list()

		#append images vertically
		#create plots
		for (col in 1:ncol(nucdf)) {

			#get peak start
			peakStart <- peakFile[col,2]
			chr <- peakFile[col,1]

			#get x axis range
			xRange <- peakStart:(peakStart + (windowSize/binsize) - 1)
			xRange <- xRange * binsize

			#modify column values
			tempNUC[1] <- xRange
			tempNFR[1] <- xRange

			tempNFR[2] <- nfrdf[col]
			tempNFR[3] <- ctrlNfrdf[col]
			tempNUC[2] <- nucdf[col]
			tempNUC[3] <- ctrlNucdf[col]

			#replace zeros to NA
			tempNFR[tempNFR == NA] <- 0
			tempNUC[tempNUC == NA] <- 0

			tempNUC.melt <- melt(tempNUC, id.vars = 'rows', variable.name = 'series')
			tempNFR.melt <- melt(tempNFR, id.vars = 'rows', variable.name = 'series')

			#define an array of type for color coding
			colorTypeNUC=rep(c("blue", "grey43"), each=(nrow(tempNUC.melt)/2))
			colorTypeNFR=rep(c("red", "grey43"), each=(nrow(tempNFR.melt)/2))

			tempNUC.melt$color <- colorTypeNUC
			tempNFR.melt$color <- colorTypeNFR

            #check if max RPM is zero in both the sample and the control. If so, ylim has to be manually set
            nucMax <- max(tempNUC.melt$value, na.rm = TRUE)
            nfrMax <- max(tempNFR.melt$value, na.rm = TRUE)

            #if max RPM in nfr bw is zero, set ylim manually
            if (nfrMax == 0) {
                NFR <- ggplot(tempNFR.melt, aes(rows,value)) +
                        geom_area(colour=tempNFR.melt$color, fill=tempNFR.melt$color) + 
                        facet_grid(series ~ .) + 
                        #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                        #labs(x = chr, y = "Score") +
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
                            #axis.title.y = element_blank()
                            #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                            #aspect.ratio=1/4
                            ) +
                        theme(strip.background = element_rect(fill="red")) +
                        theme(strip.text = element_text(colour = 'white', face="bold", size = 14)) +
                        ylim(0, 5)
                        #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            
            #if max RPM in nfr bw is not zero, let ggplot automatically set ylim
            } else {
                NFR <- ggplot(tempNFR.melt, aes(rows,value)) +
                        geom_area(colour=tempNFR.melt$color, fill=tempNFR.melt$color) + 
                        facet_grid(series ~ .) + 
                        #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                        #labs(x = chr, y = "Score") +
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
                            #axis.title.y = element_blank()
                            #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                            #aspect.ratio=1/4
                            ) +
                        theme(strip.background = element_rect(fill="red")) +
                        theme(strip.text = element_text(colour = 'white', face="bold", size = 14))
                        #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            }

            #if max RPM in nuc bw is zero, set ylim manually
            if (nucMax == 0) {
                NUC <- ggplot(tempNUC.melt, aes(rows,value)) +
                    geom_area(colour=tempNUC.melt$color, fill=tempNUC.melt$color) + 
                    facet_grid(series ~ .) + 
                    #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                    labs(x = chr, y = "Score") +
                    theme(#plot.title = element_text(hjust = 0.5, face="bold"),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA),
                        axis.title.x = element_text(face="bold"),
                        #axis.title.y = element_blank()
                        axis.title.y = element_text(colour="white")
                        #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                        #aspect.ratio=1/4
                        ) +
                    theme(strip.background = element_rect(fill="blue")) +
                    theme(strip.text = element_text(colour = 'white', face="bold", size = 14)) +
                    ylim(0, 5)
                    #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            
            #if max RPM in nuc bw is not zero:
            } else {
                NUC <- ggplot(tempNUC.melt, aes(rows,value)) +
                    geom_area(colour=tempNUC.melt$color, fill=tempNUC.melt$color) + 
                    facet_grid(series ~ .) + 
                    #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                    labs(x = chr, y = "Score") +
                    theme(#plot.title = element_text(hjust = 0.5, face="bold"),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA),
                        axis.title.x = element_text(face="bold"),
                        #axis.title.y = element_blank()
                        axis.title.y = element_text(colour="white")
                        #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                        #aspect.ratio=1/4
                        ) +
                    theme(strip.background = element_rect(fill="blue")) +
                    theme(strip.text = element_text(colour = 'white', face="bold", size = 14))
                    #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            }

            #merge both nfr and nuc visualizations vertically
			merged <- plot_grid(NFR, NUC, ncol = 1, align = "v")

            #create y-axis label
            newMerged <- ggdraw(merged) +
            draw_label("RPM", x=0.01, y=0.5, vjust=0.5, hjust=0.5, angle = 90, fontface="bold", size = 12)

            #append current sample's peak example visualization to a list
            listOfPlots <- append(listOfPlots, list(newMerged))

        } #end of loop to create plots

        #merge all peak examples
        allMerged <- plot_grid(plotlist=listOfPlots, ncol = 1, align = "v")

        #save peak examples as png and pdf
        ggsave(paste0(des, "/", "peak.examples.png"), allMerged, width = 10, height = (strtoi(numHighestPeaks) * 6), dpi = 500, units = "in", device = "png")
        ggsave(paste0(des, "/", "peak.examples.pdf"), allMerged, width = 10, height = (strtoi(numHighestPeaks) * 6), dpi = 500, units = "in", device = "pdf")

    
    #if a control sample does not exist
    #do not visualize control samples
	} else {

		#extract by 10bp
		nucPeakTable <- extractBigWigData(bed = bed, bw = nucBW, width = windowSize, binSize = binsize)
		nfrPeakTable <- extractBigWigData(bed = bed, bw = nfrBW, width = windowSize, binSize = binsize)

		#convert to DF and transpose
		nucdf <- as.data.frame(t(nucPeakTable))
		nfrdf <- as.data.frame(t(nfrPeakTable))
		tempNUC <- as.data.frame(t(nucPeakTable))
		tempNFR <- as.data.frame(t(nfrPeakTable))

		#Add x axis numbering
		tempNUC <- cbind(rows = 1:nrow(tempNUC), tempNUC)
		#remove unnecessary columns (last two)
		#this is due to errors when trying to create a new dataframe. Tedious work...
        tempNUC = tempNUC[,1:2]

		#rename column headers
		colnames(tempNUC) <- c("rows", "NUC")

		#Add x axis numbering
		tempNFR <- cbind(rows = 1:nrow(tempNFR), tempNFR)
		#remove unnecessary columns (last two)
		#this is due to errors when trying to create a new dataframe. Tedious work...
        tempNFR = tempNFR[,1:2]

		#rename column headers
		colnames(tempNFR) <- c("rows", "NFR")


		#read peak file
		peakFile <- read.table(bed, header = FALSE, sep = "\t")

        #create empty list of images
        listOfPlots <- list()

		#append images vertically
		#create plots
		for (col in 1:ncol(nucdf)) {

			#get peak start
			peakStart <- peakFile[col,2]
			chr <- peakFile[col,1]

			#get x axis range
			xRange <- peakStart:(peakStart + (windowSize/binsize) - 1)
			xRange <- xRange * binsize

			#modify column values
			tempNUC[1] <- xRange
			tempNFR[1] <- xRange

			tempNFR[2] <- nfrdf[col]
			tempNUC[2] <- nucdf[col]

			#replace zeros to NA
			tempNFR[tempNFR == NA] <- 0
			tempNUC[tempNUC == NA] <- 0

			tempNUC.melt <- melt(tempNUC, id.vars = 'rows', variable.name = 'series')
			tempNFR.melt <- melt(tempNFR, id.vars = 'rows', variable.name = 'series')

			#define an array of type for color coding
			colorTypeNUC=rep(c("blue"), each=(nrow(tempNUC.melt)))
			colorTypeNFR=rep(c("red"), each=(nrow(tempNFR.melt)))

			tempNUC.melt$color <- colorTypeNUC
			tempNFR.melt$color <- colorTypeNFR
            nucMax <- max(tempNUC.melt$value, na.rm = TRUE)
            nfrMax <- max(tempNFR.melt$value, na.rm = TRUE)

            #if max RPM in nfr bw is zero, set ylim manually
            if (nfrMax == 0) {
                NFR <- ggplot(tempNFR.melt, aes(rows,value)) +
                        geom_area(colour=tempNFR.melt$color, fill=tempNFR.melt$color) + 
                        facet_grid(series ~ .) + 
                        #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                        #labs(x = chr, y = "Score") +
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
                            #axis.title.y = element_blank()
                            #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                            #aspect.ratio=1/4
                            ) +
                        theme(strip.background = element_rect(fill="red")) +
                        theme(strip.text = element_text(colour = 'white', face="bold", size = 14)) +
                        ylim(0, 5)
                        #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            
            #if max RPM in nfr bw is zero, set ylim automatically
            } else {
                NFR <- ggplot(tempNFR.melt, aes(rows,value)) +
                        geom_area(colour=tempNFR.melt$color, fill=tempNFR.melt$color) + 
                        facet_grid(series ~ .) + 
                        #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                        #labs(x = chr, y = "Score") +
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
                            #axis.title.y = element_blank()
                            #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                            #aspect.ratio=1/4
                            ) +
                        theme(strip.background = element_rect(fill="red")) +
                        theme(strip.text = element_text(colour = 'white', face="bold", size = 14))
                        #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            }

            #if max RPM in nuc bw is zero, set ylim manually
            if (nucMax == 0) {
                NUC <- ggplot(tempNUC.melt, aes(rows,value)) +
                    geom_area(colour=tempNUC.melt$color, fill=tempNUC.melt$color) + 
                    facet_grid(series ~ .) + 
                    #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                    labs(x = chr, y = "Score") +
                    theme(#plot.title = element_text(hjust = 0.5, face="bold"),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA),
                        axis.title.x = element_text(face="bold"),
                        #axis.title.y = element_blank()
                        axis.title.y = element_text(colour="white")
                        #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                        #aspect.ratio=1/4
                        ) +
                    theme(strip.background = element_rect(fill="blue")) +
                    theme(strip.text = element_text(colour = 'white', face="bold", size = 14)) +
                    ylim(0, 5)
                    #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            
            #if max RPM in nuc bw is zero, set ylim automatically
            } else {
                NUC <- ggplot(tempNUC.melt, aes(rows,value)) +
                    geom_area(colour=tempNUC.melt$color, fill=tempNUC.melt$color) + 
                    facet_grid(series ~ .) + 
                    #labs(title=paste0("Highest Scoring Peak #", col), x = chr, y = "Score") +
                    labs(x = chr, y = "Score") +
                    theme(#plot.title = element_text(hjust = 0.5, face="bold"),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA),
                        axis.title.x = element_text(face="bold"),
                        #axis.title.y = element_blank()
                        axis.title.y = element_text(colour="white")
                        #plot.margin=grid::unit(c(0,0,0,0), "mm"),
                        #aspect.ratio=1/4
                        ) +
                    theme(strip.background = element_rect(fill="blue")) +
                    theme(strip.text = element_text(colour = 'white', face="bold", size = 14))
                    #scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
            }

            #merge nfr and nuc plots for the current sample
			merged <- plot_grid(NFR, NUC, ncol = 1, align = "v")

            #create y-axis label
            newMerged <- ggdraw(merged) +
            draw_label("RPM", x=0.01, y=0.5, vjust=0.5, hjust=0.5, angle = 90, fontface="bold", size = 12)

            #append current sample's peak example visualization to a list
            listOfPlots <- append(listOfPlots, list(newMerged))

		}

        #merge all peak examples
        allMerged <- plot_grid(plotlist=listOfPlots, ncol = 1, align = "v")

        #save peak examples as png and pdf
        ggsave(paste0(des, "/", "peak.examples.png"), allMerged, width = 10, height = (strtoi(numHighestPeaks) * 3), dpi = 500, units = "in", device = "png")
        ggsave(paste0(des, "/", "peak.examples.pdf"), allMerged, width = 10, height = (strtoi(numHighestPeaks) * 3), dpi = 500, units = "in", device = "pdf")

	}

}


#__main__
sampTemp <- head(unlist(strsplit(peakFilePath, "/")), -2)
sampDir <- paste(sampTemp, collapse = '/')

#Parse string to get sample Name
sampleName <- tail(unlist(strsplit(sampDir, "/")), n=1)

#set both nuc and nfr bw files
nucBW <- paste(sampDir, "/igv.nuc.con.bw", sep="")
nfrBW <- paste(sampDir, "/igv.nfr.con.bw", sep="")

#parse control sample
if (ctrlSampleName == "NULL") {
    ctrlNucBW <- NULL
    ctrlNfrBW <- NULL
} else {
    ctrlNucBW <- (paste0(getwd(), "/3.Sample/", ctrlSampleName, "/igv.nuc.con.bw"))
    ctrlNfrBW <- (paste0(getwd(), "/3.Sample/", ctrlSampleName, "/igv.nfr.con.bw"))
}

system(sprintf(paste0("head -n ", numHighestPeaks, " ", peakFilePath, " > ", outDir, "/topPeaks.bed")))

if (mode == 'factor') {
    windowSize = 1000
    binsize = 2
} else {
    windowSize = 10000
    binsize = 20
}

draw_peak_example(paste0(outDir, "/topPeaks.bed"), nucBW = nucBW, nfrBW = nfrBW, ctrlNucBW = ctrlNucBW, ctrlNfrBW = ctrlNfrBW, windowSize = windowSize, binsize = binsize, des=outDir)

system(sprintf(paste0("rm ", outDir, "/topPeaks.bed")))
