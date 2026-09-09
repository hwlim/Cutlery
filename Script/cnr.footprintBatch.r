#!/usr/bin/env Rscript

###############################################################################
## cnr.footprintBatch.r
##
## Stand-alone CUT&RUN / ChIP-exo style footprint analysis and visualization.
##
## Given
##   1) a BED file of regions to scan (e.g. peak calls),
##   2) a pair of 1bp-resolution stranded bigWig files (<prefix>.plus.bw / .minus.bw),
##   3) a Homer *de novo* motif search result directory (containing homerResults/),
##   4) a genome FASTA (samtools-indexed),
## this script selects enriched motifs, scans them, keeps the instances that show
## a footprint, and renders one figure per motif titled with Homer's prediction.
##
## SELF-CONTAINED BY DESIGN
##   - No source() of LimLabBase / ChipExoUtil / IDOM. All borrowed routines are
##     reimplemented below (seqToInt, drawSeqHeatmap, drawExoTemplate,
##     drawHeatmapSingle, visualizeExoAvgProfile, visualizeExoInstances,
##     extractBigWigDataStranded1bp, extractDNA, readBedFile, assertFileExist).
##   - No external binaries. Replaces:
##       bwtool extract        -> rtracklayer::import()
##       annotatePeaks.pl -m   -> Biostrings::matchPWM()
##       bedtools/homerTools   -> Biostrings::getSeq() on an indexed FASTA
##       meme2images           -> not needed (no logo panel; motif identity is
##                                carried in the figure title)
##       convert / pdf2ps      -> pdf() device rasterized by magick, the R
##                                binding to ImageMagick (no shell-out)
##   - R package deps are CRAN/Bioconductor only, declared and checked up front:
##       optparse, magick, Biostrings, rtracklayer, GenomicRanges, Rsamtools,
##       parallel  (parallel ships with R; gplots/fields/PWMEnrich are NOT used)
##
## NOTE: magick reads PDF through ImageMagick's Ghostscript delegate. If PDF
## rasterization fails, the ImageMagick build lacks that delegate.
##
## HOMER SCORING CONVENTION
##   Homer .motif files hold a probability matrix (rows = positions, cols = ACGT)
##   and a log-odds detection threshold in header field 3. That threshold is on a
##   NATURAL LOG scale against a uniform 0.25 background. Verified empirically over
##   1313 motifs from the Homer database: under natural log the threshold is always
##   reachable (median 0.63 of max score); under log10 it is unreachable for 96% of
##   motifs. Biostrings works in log2, so the matrix is converted with log2() and
##   the threshold is divided by log(2).
##
##   NOTE: this reproduces Homer's documented convention but has not been
##   validated hit-for-hit against `annotatePeaks.pl -mbed`. If exact agreement
##   matters, run Homer once on a fixed input and compare (see --minScoreScale).
##
##      Written for the Lim Lab. Derived from cnr.analyzeFootprintBatch.r
##      (Hee Woong Lim) and idom.visualizeExoBed.r.
###############################################################################

###############################################################################
## Required packages
##
## Declared up front and checked before anything else runs, so a missing
## dependency reports itself immediately with an install hint rather than
## failing part way through an analysis.
###############################################################################

requiredPackages <- c(
	optparse      = "CRAN",
	magick        = "CRAN",
	Biostrings    = "Bioconductor",
	rtracklayer   = "Bioconductor",
	GenomicRanges = "Bioconductor",
	Rsamtools     = "Bioconductor",
	seqLogo       = "Bioconductor",
	grid          = "base",
	parallel      = "base"
)

local({
	have <- vapply(names(requiredPackages),
				   function(p) requireNamespace(p, quietly = TRUE), logical(1))
	if (all(have)) return(invisible(NULL))

	missing <- names(requiredPackages)[!have]
	write("Error: required R package(s) are not available:", stderr())
	for (p in missing)
		write(sprintf("  - %-14s (%s)", p, requiredPackages[[p]]), stderr())

	cran <- missing[requiredPackages[missing] == "CRAN"]
	bioc <- missing[requiredPackages[missing] == "Bioconductor"]
	write("", stderr())
	write("Install with:", stderr())
	if (length(cran) > 0)
		write(sprintf("  install.packages(c(%s))",
					  paste(sprintf('"%s"', cran), collapse = ", ")), stderr())
	if (length(bioc) > 0) {
		write('  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")', stderr())
		write(sprintf("  BiocManager::install(c(%s))",
					  paste(sprintf('"%s"', bioc), collapse = ", ")), stderr())
	}
	if ("parallel" %in% missing)
		write("  'parallel' ships with R; this R installation looks incomplete.", stderr())
	quit(save = "no", status = 1)
})

suppressPackageStartupMessages({
	library(optparse)
	library(magick)
	library(Biostrings)
	library(rtracklayer)
	library(GenomicRanges)
	library(Rsamtools)
	library(seqLogo)
	library(grid)
	library(parallel)
})


###############################################################################
## Command line
###############################################################################

option_list <- list(
	make_option(c("-o","--outPrefix"), default="footprint",
		help="Output prefix (may include a path). default=footprint"),
	make_option(c("-n","--name"), default=NULL,
		help="Data set name used in figure titles. default=<basename of bed>"),
	make_option(c("-g","--genome"), default=NULL,
		help="(required) Genome FASTA. A samtools .fai index must exist next to it"),
	make_option(c("-t","--minTargetPercent"), default=10,
		help="Minimum Homer target %% for motif selection. default=10"),
	make_option(c("-P","--maxPvalue"), default=1e-10,
		help="Maximum Homer motif p-value for selection. default=1e-10"),
	make_option(c("-m","--maxMotifCount"), default=10,
		help="Maximum number of Homer motifs to consider. default=10"),
	make_option(c("--margin"), default=20,
		help="Flank margin (bp) used for the footprint contrast test. default=20"),
	make_option(c("--seqMargin"), default=20,
		help="Margin (bp) around the motif for the DNA sequence panel. default=20"),
	make_option(c("--exoMargin"), default=50,
		help="Margin (bp) around the motif for the signal panels. default=50"),
	make_option(c("--minContrast"), default=0,
		help="Keep motif instances with footprint contrast strictly above this. default=0"),
	make_option(c("--maxInstances"), default=0,
		help="Cap instances drawn per motif, keeping the highest contrast. 0 = no cap. default=0"),
	make_option(c("--ymax"), default="0,0",
		help="Comma-separated y-max for average profile and heatmap. 0 = automatic. default=0,0"),
	make_option(c("-c","--combineMode"), default="sum",
		help="How to combine plus/minus for the heatmap: sum or max. default=sum"),
	make_option(c("--minScoreScale"), default=1.0,
		help="Multiplier applied to the converted Homer threshold. <1 loosens, >1 tightens. default=1.0"),
	make_option(c("--density"), default=200,
		help="Raster density (DPI) used when converting the PDF to PNG. default=200"),
	make_option(c("--pdf"), default=FALSE, action="store_true",
		help="Keep the intermediate PDF next to the PNG. default=FALSE"),
	make_option(c("--allowEmpty"), default=FALSE, action="store_true",
		help="Treat 'nothing to show' as success: if no motif can be read, passes selection,
		is found in the regions, or yields a footprint, write <outPrefix>.4.complete and exit 0
		instead of failing. Genuine errors (missing files, a motif that errors during
		processing) still fail. Intended for pipeline use. default=FALSE"),
	make_option(c("-p","--parallel"), default=1,
		help="Number of cores for per-motif processing. default=1"),
	make_option(c("-v","--verbose"), default=FALSE, action="store_true",
		help="Verbose progress messages")
)

parser <- OptionParser(
	usage = "%prog [options] -g <genome.fa> <regions.bed> <homer motif dir> <bigWig prefix>",
	option_list = option_list,
	description = "Description:
	Footprint analysis and visualization for 1bp-resolution stranded coverage.

Input:
	<regions.bed>       BED file of regions to scan for motifs (e.g. peak calls)
	<homer motif dir>   Homer de novo result dir; must contain homerResults/motif<i>.motif
	<bigWig prefix>     Prefix such that <prefix>.plus.bw and <prefix>.minus.bw exist
	-g <genome.fa>      Indexed genome FASTA used for sequence extraction and scanning

Output:
	<outPrefix>.0.selectedMotif.txt      selected motifs and their Homer statistics
	<outPrefix>.1.motifScan.bed          all motif instances found in the regions
	<outPrefix>.3.motif<NN>.<motif>/     one folder per motif, Homer rank first
	    CnR.viz.png                      sequence + average profile + signal heatmap
	    CnR.viz.avg.png                  average profile alone
	    CnR.sorted.bed                   instances passing the contrast test,
	                                     in plotted order, with Contrast/Anchor/Flank
	    CnR.sorted.fa                    sequences in plotted order
	<outPrefix>.5.summary.txt            per-motif status and instance counts
	<outPrefix>.4.complete               written only if every motif succeeded")

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args) != 3) {
	print_help(parser)
	stop("Requires exactly three positional arguments: <regions.bed> <homer motif dir> <bigWig prefix>", call. = FALSE)
}
src.bed      <- arguments$args[1]
src.motifDir <- arguments$args[2]
src.bwPrefix <- arguments$args[3]

if (is.null(opt$genome)) stop("Genome FASTA (-g/--genome) is required", call. = FALSE)
if (!opt$combineMode %in% c("sum", "max")) stop("--combineMode must be 'sum' or 'max'", call. = FALSE)

outPrefix <- opt$outPrefix
dataName  <- if (is.null(opt$name)) sub("\\.bed(\\.gz)?$", "", basename(src.bed)) else opt$name

ymaxL <- as.numeric(strsplit(opt$ymax, ",")[[1]])
if (length(ymaxL) != 2 || any(is.na(ymaxL)) || any(ymaxL < 0))
	stop("--ymax must be two non-negative comma-separated numbers", call. = FALSE)
yMaxAvg     <- ymaxL[1]
yMaxHeatmap <- ymaxL[2]

src.bwPlus  <- sprintf("%s.plus.bw",  src.bwPrefix)
src.bwMinus <- sprintf("%s.minus.bw", src.bwPrefix)


###############################################################################
## Small utilities (internalized from LimLabBase)
###############################################################################

assertFileExist <- function(paths) {
	missing <- paths[!file.exists(paths)]
	if (length(missing) > 0)
		stop(sprintf("File does not exist: %s", paste(missing, collapse = ", ")), call. = FALSE)
	invisible(TRUE)
}

vmsg <- function(fmt, ...) {
	if (isTRUE(opt$verbose)) write(sprintf(fmt, ...), stderr())
}

msg <- function(fmt, ...) write(sprintf(fmt, ...), stderr())

## Read a BED file into a data.frame with at least 6 columns.
readBed <- function(path) {
	d <- read.delim(path, header = FALSE, stringsAsFactors = FALSE,
					comment.char = "#", blank.lines.skip = TRUE)
	if (ncol(d) < 3) stop(sprintf("%s does not look like a BED file", path), call. = FALSE)
	if (ncol(d) < 4) d[[4]] <- sprintf("R%d", seq_len(nrow(d)))
	if (ncol(d) < 5) d[[5]] <- 0
	if (ncol(d) < 6) d[[6]] <- "+"
	d <- d[, 1:6]
	colnames(d) <- c("chr", "start", "end", "name", "score", "strand")
	d$chr    <- as.character(d$chr)
	d$name   <- as.character(d$name)
	d$strand <- as.character(d$strand)
	d$strand[!d$strand %in% c("+", "-")] <- "+"
	d
}

writeBed <- function(df, path) {
	write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

## BED (0-based, half-open) -> GRanges (1-based, inclusive)
bedToGRanges <- function(d) {
	GRanges(seqnames = d$chr,
			ranges   = IRanges(start = d$start + 1, end = d$end),
			strand   = d$strand,
			name     = d$name)
}


###############################################################################
## Homer de novo motif parsing
###############################################################################

## Parse the native header line of a Homer .motif file, e.g.
##   >AWRACAAWRG  1-AWRACAAWRG,BestGuess:SOX15/MA1152.1/Jaspar(0.974)  7.5  -1234.5  0  T:100.0(20.00%),B:50.0(2.00%),P:1e-50
## Handles the awkward name cases: "Pou5f1::Sox2" and "TEAD(TEA)/Fibroblast-..."
parseHomerMotifHeader <- function(hdr) {
	el <- strsplit(hdr, "\t")[[1]]
	res <- list(consensus = sub("^>", "", el[1]), raw = hdr)

	## Best-guess name
	nm <- NA_character_
	if (length(el) >= 2) {
		part <- strsplit(el[2], ",BestGuess:", fixed = TRUE)[[1]]
		if (length(part) >= 2) {
			nm <- strsplit(part[2], "[/()]")[[1]][1]
		} else {
			nm <- el[2]
		}
		nm <- gsub("::", "_", nm)
	}
	if (is.na(nm) || !nzchar(nm)) nm <- res$consensus
	## Make the name safe for use as a file / directory component
	res$name <- gsub("[^A-Za-z0-9_.+-]", "_", nm)

	res$score <- suppressWarnings(as.numeric(el[3]))
	res$logp  <- suppressWarnings(as.numeric(el[4]))

	## Statistics field: T:####(##%),B:####(##%),P:1e-##
	stat <- if (length(el) >= 6) el[6] else ""
	grab <- function(pat, s) {
		m <- regmatches(s, regexpr(pat, s, perl = TRUE))
		if (length(m) == 0) return(NA_real_)
		suppressWarnings(as.numeric(sub(pat, "\\1", m, perl = TRUE)))
	}
	res$target     <- grab("T:[0-9.eE+-]+\\(([0-9.]+)%\\)", stat)
	res$background <- grab("B:[0-9.eE+-]+\\(([0-9.]+)%\\)", stat)
	res$pvalue     <- grab("P:([0-9.eE+-]+)", stat)
	res
}

## Read homerResults/motif<i>.motif for i in 1..maxCount
##
## A missing motifDir is a path error and always fails. A motifDir that exists
## but has no homerResults/ is what Homer leaves behind when it found nothing
## (run_homer_motif touches homerResults.html without creating the subdirectory),
## so that returns no motifs and is left for the caller to handle.
readHomerMotifDir <- function(motifDir, maxCount) {
	if (!dir.exists(motifDir))
		stop(sprintf("Motif directory does not exist: %s", motifDir), call. = FALSE)
	resultDir <- file.path(motifDir, "homerResults")
	if (!dir.exists(resultDir)) {
		msg("Warning: %s has no homerResults/ -- treating as no motifs found", motifDir)
		return(list())
	}

	out <- list()
	for (i in seq_len(maxCount)) {
		f <- file.path(resultDir, sprintf("motif%d.motif", i))
		if (!file.exists(f)) break
		lines <- readLines(f, warn = FALSE)
		lines <- lines[nzchar(lines)]
		hdrIdx <- grep("^>", lines)
		if (length(hdrIdx) == 0) {
			msg("  Warning: no header in %s, skipping", f)
			next
		}
		info <- parseHomerMotifHeader(lines[hdrIdx[1]])
		## matrix rows run until the next header (or end of file)
		last <- if (length(hdrIdx) > 1) hdrIdx[2] - 1 else length(lines)
		body <- lines[(hdrIdx[1] + 1):last]
		mat  <- do.call(rbind, lapply(strsplit(body, "[\t ]+"), as.numeric))
		if (is.null(dim(mat)) || ncol(mat) != 4) {
			msg("  Warning: unexpected matrix shape in %s, skipping", f)
			next
		}
		colnames(mat) <- c("A", "C", "G", "T")
		info$rank <- i
		info$prob <- mat
		info$file <- f
		out[[length(out) + 1]] <- info
	}
	out
}

## Homer probability matrix -> Biostrings-style log2 odds PWM (4 x L, rownames ACGT)
homerToLog2Pwm <- function(prob, bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25), pseudo = 1e-3) {
	p <- pmax(prob, pseudo)
	p <- p / rowSums(p)
	lo <- log2(sweep(p, 2, bg[colnames(p)], "/"))
	pwm <- t(lo)
	rownames(pwm) <- colnames(p)
	pwm <- pwm[c("A", "C", "G", "T"), , drop = FALSE]
	storage.mode(pwm) <- "double"
	pwm
}


###############################################################################
## Sequence extraction (replaces extractDNA / bedtools / homerTools)
###############################################################################

extractSeq <- function(gr, faFile) {
	## getSeq on FaFile ignores strand, so revcomp manually where needed
	gr0 <- gr
	strand(gr0) <- "*"
	s <- getSeq(faFile, gr0)
	neg <- which(as.character(strand(gr)) == "-")
	if (length(neg) > 0) s[neg] <- reverseComplement(s[neg])
	names(s) <- mcols(gr)$name
	s
}


###############################################################################
## BigWig extraction (replaces extractBigWigDataStranded1bp / bwtool)
###############################################################################

## Returns list(plus=matrix, minus=matrix): rows = regions (order of gr), cols = bp.
## Positions with no data in the bigWig come back as 0, matching `bwtool -fill=0`.
## flipByStrand mirrors readExoProfile2(): for minus-strand regions the two
## strands are swapped, reversed, and negated so that every profile is oriented
## with respect to the motif.
extractStranded1bp <- function(gr, bwPlus, bwMinus, flipByStrand = TRUE) {
	w <- unique(width(gr))
	if (length(w) != 1) stop("extractStranded1bp requires equal-width regions", call. = FALSE)

	gr0 <- gr
	strand(gr0) <- "*"

	readOne <- function(path) {
		bwf <- BigWigFile(path)
		## A range on a contig absent from the bigWig is a hard error in
		## rtracklayer (bwtool would have zero-filled), so guard explicitly.
		keep <- as.character(seqnames(gr0)) %in% seqlevels(seqinfo(bwf))
		m <- matrix(0, nrow = length(gr0), ncol = w)
		if (any(keep)) {
			nl <- import(bwf, which = gr0[keep], as = "NumericList")
			## NOTE: nl[, idx] is a silent no-op on NumericList; always go via a matrix.
			m[keep, ] <- matrix(as.numeric(unlist(nl)), nrow = sum(keep), byrow = TRUE)
		}
		m
	}

	mP <- readOne(bwPlus)
	mM <- readOne(bwMinus)

	if (flipByStrand) {
		neg <- which(as.character(strand(gr)) == "-")
		if (length(neg) > 0) {
			tmpP <- mP[neg, , drop = FALSE]
			mP[neg, ] <- -mM[neg, ncol(mM):1, drop = FALSE]
			mM[neg, ] <- -tmpP[, ncol(tmpP):1, drop = FALSE]
		}
	}
	list(plus = mP, minus = mM)
}


###############################################################################
## Motif scanning (replaces annotatePeaks.pl -m -mbed)
###############################################################################

## Scan a PWM against region sequences on both strands and return a GRanges of
## genomic motif instances. Sequences are concatenated with runs of N so a single
## matchPWM() call covers every region; N never matches, so no hit can straddle
## a boundary (asserted below).
scanMotif <- function(seqs, gr, pwm, minScore) {
	L <- ncol(pwm)
	if (any(width(gr) < L)) {
		keep <- width(gr) >= L
		seqs <- seqs[keep]; gr <- gr[keep]
	}
	if (length(seqs) == 0) return(GRanges())

	spacer   <- paste(rep("N", L), collapse = "")
	widths   <- width(seqs)
	offsets  <- cumsum(c(0, head(widths + L, -1))) + 1   # 1-based start of each region
	big      <- DNAString(paste(paste0(as.character(seqs), spacer), collapse = ""))

	collect <- function(p, strandChar) {
		## matchPWM warns that the N spacers are not in [ACGT] and gives them weight 0.
		## That is expected: the spacers exist only to separate regions, and any hit
		## that reaches into one is discarded by the boundary check below.
		v <- withCallingHandlers(
			matchPWM(p, big, min.score = minScore),
			warning = function(w) {
				if (grepl("letters not in \\[ACGT\\]", conditionMessage(w))) invokeRestart("muffleWarning")
			})
		if (length(v) == 0) return(NULL)
		st <- start(v)
		idx <- findInterval(st, offsets)
		within <- st - offsets[idx] + 1
		ok <- within >= 1 & (within + L - 1) <= widths[idx]
		if (!all(ok)) {   # should be impossible; N blocks cross-boundary hits
			idx <- idx[ok]; within <- within[ok]
		}
		if (length(idx) == 0) return(NULL)
		gStart <- start(gr)[idx] + within - 1
		GRanges(seqnames = seqnames(gr)[idx],
				ranges   = IRanges(start = gStart, width = L),
				strand   = strandChar,
				region   = mcols(gr)$name[idx])
	}

	hits <- c(collect(pwm, "+"), collect(reverseComplement(pwm), "-"))
	if (is.null(hits) || length(hits) == 0) return(GRanges())
	sort(hits)
}


###############################################################################
## Footprint contrast
###############################################################################

## log2( flank signal / motif signal ). Positive means the motif is protected
## relative to its flanks, i.e. a footprint. Strand is deliberately ignored here
## (the window is symmetric), matching the original implementation.
footprintContrast <- function(gr, bwPlus, bwMinus, margin = 20, pseudo = 0.1) {
	w <- unique(width(gr))
	if (length(w) != 1) stop("footprintContrast requires equal-width regions", call. = FALSE)

	ext <- gr
	ranges(ext) <- IRanges(start = start(gr) - margin, end = end(gr) + margin)
	prof <- extractStranded1bp(ext, bwPlus, bwMinus, flipByStrand = FALSE)

	ia <- (margin + 1):(margin + w)
	fl <- c(seq_len(margin), (margin + w + 1):(2 * margin + w))

	sigAnchor <- rowMeans(prof$plus[, ia, drop = FALSE]) + abs(rowMeans(prof$minus[, ia, drop = FALSE]))
	sigFlank  <- rowMeans(prof$plus[, fl, drop = FALSE]) + abs(rowMeans(prof$minus[, fl, drop = FALSE]))

	data.frame(Anchor   = sigAnchor,
			   Flank    = sigFlank,
			   Contrast = log2((sigFlank + pseudo) / (sigAnchor + pseudo)))
}


###############################################################################
## Drawing (internalized from genomeR.r / commonR.r / ChipExoUtil.r)
###############################################################################

## Render a figure by drawing into a PDF device and rasterizing it with magick
## (the R binding to ImageMagick). Drawing as vector first and rasterizing at a
## chosen density gives sharper text and thinner lines than a direct png()
## device, and reproduces the original `pdf()` + `convert -density 200` route
## without shelling out.
##
## width/height are in INCHES (pdf() units). `expr` is evaluated lazily, after
## the device is open.
renderFigure <- function(outPng, width, height, expr, density = 200, keepPdf = FALSE) {
	pdfPath <- sub("\\.png$", ".pdf", outPng)
	pdf(pdfPath, width = width, height = height, bg = "white")
	tryCatch(expr, finally = invisible(dev.off()))

	## image_read() rasterizes the PDF through ImageMagick's own Ghostscript
	## delegate -- the same path the original `convert -density 200` used.
	## (image_read_pdf() is deliberately avoided: it routes through the extra
	## pdftools/poppler dependency instead.)
	img <- tryCatch(
		magick::image_read(pdfPath, density = density),
		error = function(e) stop(sprintf(
			"Could not rasterize %s (%s).\n  ImageMagick needs its Ghostscript delegate to read PDF; check that `gs` is installed.",
			pdfPath, conditionMessage(e)), call. = FALSE))
	## Take the last page: grid-based drawing (seqLogo) calls grid.newpage(),
	## which can leave an empty leading page. Base-graphics figures are
	## single-page, so last == only.
	if (length(img) > 1) img <- img[length(img)]
	magick::image_write(img, path = outPng, format = "png")

	if (!keepPdf) unlink(pdfPath)
	invisible(outPng)
}

## "ACGT" -> 0,1,2,3 as an integer matrix, one row per sequence.
seqToInt <- function(seqs) {
	ch <- do.call(rbind, strsplit(as.character(seqs), "", fixed = TRUE))
	out <- matrix(NA_integer_, nrow(ch), ncol(ch))
	out[ch == "A" | ch == "a"] <- 0L
	out[ch == "C" | ch == "c"] <- 1L
	out[ch == "G" | ch == "g"] <- 2L
	out[ch == "T" | ch == "t"] <- 3L
	out
}

drawSeqHeatmap <- function(xtick, seqs, main = "", mar = c(3, 2, 2, 2), cex = 1) {
	## A red, C blue, G yellow, T green (matching the original palette)
	pal <- c("red", "blue", "yellow", "green")
	sm  <- seqToInt(seqs)
	par(mar = mar)
	image(xtick, y = seq_len(nrow(sm)), t(sm[nrow(sm):1, , drop = FALSE]),
		  col = pal, axes = FALSE, xlab = "", ylab = "", main = main)
	box()
	axis(1, cex.axis = cex)
}

## Information-content sequence logo, drawn from the DNA of the retained motif
## instances rather than from Homer's PWM -- so it shows what was actually found
## in this dataset.
##
## Colours follow the MEME convention used by drawSeqHeatmap() (A red, C blue,
## G yellow, T green), with two deliberate substitutions for legibility as
## letterforms on a white page, where the heatmap's dense adjacent cells are not
## there to carry the colour:
##   G  gold  #F5C710  instead of pure yellow
##   T  darkgreen      instead of green
## The heatmap palette in drawSeqHeatmap() is intentionally left unchanged.
logoFill <- c(A = "red", C = "blue", G = "#F5C710", T = "darkgreen")

## Height/width of the combined figure's upper-left panel: the figure is 6x8 in
## with layout heights c(1.5, 5), so that panel is about 3.0 x 1.85 in. The embed
## is rendered near this aspect so it fills the panel rather than letterboxing.
logoPanelAspect <- 0.62

## Fraction of the panel the embedded logo occupies. 1.0 fills it edge to edge,
## which crowds the neighbouring panels; lower this to shrink the logo.
logoPanelScale <- 0.78

## Note: seqLogo() calls grid.newpage() itself and labels the x-axis 1..L, so the
## title is added afterwards (it pops back to the root viewport, leaving a 2-line
## top margin) and the window is described in the title rather than the axis.
## fontsize: seqLogo derives its margins from this (2 + size/3.5 lines) and
## places the axis titles 3 lines out, so the 15pt default overflows and clips
## the y-axis label on a figure this size.
drawSeqLogo <- function(seqs, title = "", fontsize = 9) {
	cm <- consensusMatrix(seqs, as.prob = TRUE, baseOnly = TRUE)
	cm <- cm[c("A", "C", "G", "T"), , drop = FALSE]
	cs <- colSums(cm)
	cs[cs == 0] <- 1                     # columns that are entirely non-ACGT
	cm <- sweep(cm, 2, cs, "/")
	seqLogo::seqLogo(seqLogo::makePWM(cm), ic.scale = TRUE, fill = logoFill,
					 xfontsize = fontsize, yfontsize = fontsize)
	if (nzchar(title))
		grid::grid.text(title,
						y  = grid::unit(1, "npc") - grid::unit(0.8, "lines"),
						gp = grid::gpar(fontsize = 9, fontface = "bold"))
	invisible(cm)
}

## blue -> white -> red, replacing gplots::bluered
blueredPal <- colorRampPalette(c("blue", "white", "red"))

drawHeatmapSingle <- function(data, breaks = NULL, colorFun = blueredPal,
							  drawBox = FALSE, main = "", margin = c(5.1, 4.1, 4.1, 2.1)) {
	data <- as.matrix(data)
	data <- t(data[nrow(data):1, , drop = FALSE])
	if (is.null(breaks)) {
		maxX <- max(data); minX <- min(data)
		if (!is.finite(maxX) || !is.finite(minX) || maxX == minX) {
			maxX <- minX + 1
		}
		breaks <- seq(minX, maxX, by = (maxX - minX) / 100)
	}
	maxVal <- max(breaks); minVal <- min(breaks)
	data[data > maxVal] <- maxVal
	data[data < minVal] <- minVal
	par(mar = margin)
	image(data, col = colorFun(length(breaks) - 1), breaks = breaks, axes = FALSE, main = main)
	if (drawBox) box()
}

## Average plus/minus profile with the motif footprint marked at the center.
drawExoTemplate <- function(exo, yMax = 0, motifLen = 0) {
	width <- nrow(exo)
	if (motifLen > 0) {
		if ((width - motifLen) %% 2 != 0)
			stop("Profile width and motif length have different parity", call. = FALSE)
		mhalf <- (motifLen - 1) / 2
	}
	if (yMax == 0) yMax <- max(abs(exo)) * 1.05
	if (!is.finite(yMax) || yMax <= 0) yMax <- 1

	half <- (width - 1) / 2
	x <- seq(-half, half, by = 1)

	par(las = 1)
	matplot(x, exo, type = "l", lty = 1, ylim = c(-yMax, yMax), xlim = c(min(x), max(x)),
			xaxs = "i", yaxs = "i", col = c("firebrick1", "royalblue"),
			xlab = "Distance from Motif Center (bp)", ylab = "Average Signal (RPM)")
	polygon(c(x, rev(x)), c(exo[, "plus"],  rep(0, width)), col = "firebrick1", border = NA)
	polygon(c(x, rev(x)), c(exo[, "minus"], rep(0, width)), col = "royalblue",  border = NA)
	box()

	if (motifLen > 0) {
		my <- yMax / 30
		polygon(c(-mhalf - 0.5, mhalf + 0.5, mhalf + 0.5, -mhalf - 0.5),
				c(my, my, -my, -my), border = NA, col = "black")
	}
}

visualizeAvgProfile <- function(fwd, rev_, title = "", motifLen = 0, yMaxAvg = 0) {
	mat <- cbind(plus = colMeans(fwd), minus = colMeans(rev_))
	drawExoTemplate(mat, yMax = yMaxAvg, motifLen = motifLen)
	title(title)
}

## Three-panel figure: DNA sequence heatmap, average profile, signal heatmap.
## logoRaster: optional raster of the motif logo, drawn into the otherwise empty
## upper-left panel. It arrives pre-rendered because seqLogo draws through grid,
## which cannot target a base layout() panel; embedding it as a raster avoids
## pulling in gridBase.
visualizeInstances <- function(seqs, fwd, rev_, title = "", motifLen = 0,
							   yMaxAvg = 0, yMaxHeatmap = 0, cex = 1, mode = "sum",
							   logoRaster = NULL) {
	if (yMaxHeatmap == 0) {
		tmp <- apply(cbind(fwd, -rev_), 1, max)
		yMaxHeatmap <- stats::median(tmp) * 1.1
	}
	if (!is.finite(yMaxHeatmap) || yMaxHeatmap <= 0) yMaxHeatmap <- 1

	par(las = 1, oma = c(0, 0, 5, 0), cex = cex)
	layout(matrix(1:4, 2, 2, byrow = FALSE), widths = c(1, 1), heights = c(1.5, 5))

	## Panel 1 (upper left): motif logo, scaled to fit without distortion.
	## Margins kept minimal -- the raster is pre-trimmed, so this panel needs
	## only enough breathing room to keep the logo off the panel edges.
	par(mar = c(0.3, 0.3, 0.3, 0.3))
	plot.new()
	if (!is.null(logoRaster)) {
		pin <- par("pin")
		imgAsp <- nrow(logoRaster) / ncol(logoRaster)
		panAsp <- pin[2] / pin[1]
		if (imgAsp > panAsp) { h <- 1; w <- panAsp / imgAsp } else { w <- 1; h <- imgAsp / panAsp }
		w <- w * logoPanelScale
		h <- h * logoPanelScale
		rasterImage(logoRaster, 0.5 - w / 2, 0.5 - h / 2, 0.5 + w / 2, 0.5 + h / 2,
					interpolate = TRUE)
	}

	seqMargin <- (nchar(as.character(seqs)[1]) - motifLen) / 2
	exoMargin <- (ncol(fwd) - motifLen) / 2
	mtext(sprintf("%s\nyMaxHeatmap = %.3f / Margin = %d, %d (bp)",
				  title, yMaxHeatmap, seqMargin, exoMargin), outer = TRUE, line = 0.5)

	## 1) consensus sequence heatmap
	w <- nchar(as.character(seqs)[1])
	xtick.seq <- seq(-(w - 1) / 2, (w - 1) / 2, by = 1)
	drawSeqHeatmap(xtick.seq, seqs, mar = c(3, 2, 2, 2))

	## 2) average profile
	mat <- cbind(plus = colMeans(fwd), minus = colMeans(rev_))
	par(mar = c(1, 2, 2, 2))
	drawExoTemplate(mat, yMax = yMaxAvg, motifLen = motifLen)

	## 3) per-instance signal heatmap
	by <- yMaxHeatmap / 20
	if (mode == "sum") {
		allSig <- fwd + rev_
	} else {
		allSig <- ifelse(fwd > -rev_, fwd, rev_)
	}
	drawHeatmapSingle(allSig, breaks = seq(-yMaxHeatmap, yMaxHeatmap, by),
					  margin = c(3, 2, 2, 2), drawBox = TRUE)
}


###############################################################################
## Main
###############################################################################

assertFileExist(c(src.bed, src.bwPlus, src.bwMinus, opt$genome))
if (!file.exists(paste0(opt$genome, ".fai")))
	stop(sprintf("Missing FASTA index: %s.fai  (run: samtools faidx %s)", opt$genome, opt$genome), call. = FALSE)

outDir <- dirname(outPrefix)
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

## "Nothing to visualize" outcome. Under --allowEmpty this is a successful run
## with no output: write an empty summary plus the completion flag and exit 0, so
## a weak sample cannot break a pipeline DAG. Otherwise it is an error.
## Genuine errors never route through here.
finishEmpty <- function(reason) {
	if (!isTRUE(opt$allowEmpty)) stop(reason, call. = FALSE)
	msg("%s", reason)
	msg("--allowEmpty set: writing completion flag and exiting successfully")
	write.table(data.frame(Motif = character(0), Status = character(0), Instances = integer(0)),
				sprintf("%s.5.summary.txt", outPrefix),
				sep = "\t", quote = FALSE, row.names = FALSE)
	writeLines("", sprintf("%s.4.complete", outPrefix))
	quit(save = "no", status = 0)
}

msg("=================================")
msg("CUT&RUN footprint analysis")
msg("=================================")
msg("Name          = %s", dataName)
msg("Regions       = %s", src.bed)
msg("Homer motifs  = %s", src.motifDir)
msg("bigWig prefix = %s", src.bwPrefix)
msg("Genome        = %s", opt$genome)
msg("")

fa <- FaFile(opt$genome)
open(fa)
on.exit(try(close(fa), silent = TRUE), add = TRUE)
genomeSeqlevels <- as.character(seqnames(scanFaIndex(fa)))

## ---------------------------------------------------------------- regions ---
bed <- readBed(src.bed)
if (any(duplicated(bed$name))) bed$name <- sprintf("R%d", seq_len(nrow(bed)))
regions <- bedToGRanges(bed)
regions <- regions[as.character(seqnames(regions)) %in% genomeSeqlevels]
if (length(regions) == 0) stop("No regions remain after matching against the genome index", call. = FALSE)
msg("Regions to scan: %d", length(regions))

## ------------------------------------------------------------ motif input ---
motifs <- readHomerMotifDir(src.motifDir, opt$maxMotifCount)
if (length(motifs) == 0) finishEmpty("No Homer motifs could be read")

selected <- list()
for (m in motifs) {
	tag <- sprintf("motif%d (%s)", m$rank, m$name)
	if (!is.na(m$pvalue) && m$pvalue > opt$maxPvalue) {
		msg("  skip %-28s p-value %.3g > %.3g", tag, m$pvalue, opt$maxPvalue); next
	}
	if (!is.na(m$target) && m$target < opt$minTargetPercent) {
		msg("  skip %-28s target %.2f%% < %g%%", tag, m$target, opt$minTargetPercent); next
	}
	if (is.na(m$score)) {
		msg("  skip %-28s no log-odds threshold in header", tag); next
	}
	selected[[length(selected) + 1]] <- m
	msg("  keep %-28s len=%2d  thr=%.3f  p=%.3g  T=%.2f%%",
		tag, nrow(m$prob), m$score, m$pvalue, m$target)
}
if (length(selected) == 0) finishEmpty("No motif passed the selection criteria")

## Disambiguate names so output paths never collide
nameCount <- table(vapply(selected, function(x) x$name, ""))
for (i in seq_along(selected)) {
	if (nameCount[[selected[[i]]$name]] > 1)
		selected[[i]]$name <- sprintf("%s.m%d", selected[[i]]$name, selected[[i]]$rank)
}

sel.tab <- do.call(rbind, lapply(selected, function(m) data.frame(
	Rank = m$rank, Name = m$name, Consensus = m$consensus, Length = nrow(m$prob),
	LogOddsThreshold = m$score, LogP = m$logp, Pvalue = m$pvalue,
	TargetPct = m$target, BackgroundPct = m$background, stringsAsFactors = FALSE)))
write.table(sel.tab, sprintf("%s.0.selectedMotif.txt", outPrefix),
			sep = "\t", quote = FALSE, row.names = FALSE)
msg("")

## ------------------------------------------------------- region sequences ---
msg("1) Extracting region sequences")
regionSeq <- extractSeq(regions, fa)

## ------------------------------------------------------------------- scan ---
msg("2) Scanning %d motif(s)", length(selected))
scanL <- list()
for (m in selected) {
	pwm      <- homerToLog2Pwm(m$prob)
	minScore <- (m$score / log(2)) * opt$minScoreScale   # Homer natural log -> log2
	hits     <- scanMotif(regionSeq, regions, pwm, minScore)
	if (length(hits) > 0) mcols(hits)$motif <- m$name
	scanL[[m$name]] <- hits
	msg("   %-28s min.score=%7.3f (log2)   hits=%d", m$name, minScore, length(hits))
}

allHits <- suppressWarnings(do.call(c, unname(scanL[lengths(scanL) > 0])))
if (is.null(allHits) || length(allHits) == 0) finishEmpty("No motif instances found in the given regions")
writeBed(data.frame(chr = as.character(seqnames(allHits)),
					start = start(allHits) - 1, end = end(allHits),
					name = mcols(allHits)$motif, score = 0,
					strand = as.character(strand(allHits))),
		 sprintf("%s.1.motifScan.bed", outPrefix))
msg("")

## ------------------------------------------------------- per-motif figure ---
msg("3) Footprint contrast and visualization")

processMotif <- function(m) {
	hits <- scanL[[m$name]]
	motifLen <- nrow(m$prob)
	label <- m$name
	## Own FASTA handle: a handle opened in the parent is not safe to reuse
	## across a fork() from mclapply().
	faLocal <- FaFile(opt$genome)
	if (is.null(hits) || length(hits) == 0) {
		msg("   %-28s no instances, skipped", label)
		return(list(name = label, status = "no_instances", n = 0))
	}

	## Drop instances whose extended windows would fall off the contig
	pad <- max(opt$margin, opt$seqMargin, opt$exoMargin)
	si  <- seqinfo(faLocal)
	chrLen <- setNames(as.numeric(seqlengths(si)), seqlevels(si))
	ok <- (start(hits) - pad) >= 1 & (end(hits) + pad) <= chrLen[as.character(seqnames(hits))]
	ok[is.na(ok)] <- FALSE
	hits <- hits[ok]
	if (length(hits) == 0) {
		msg("   %-28s all instances too close to contig ends, skipped", label)
		return(list(name = label, status = "no_instances", n = 0))
	}

	## Footprint contrast test
	ct <- footprintContrast(hits, src.bwPlus, src.bwMinus, margin = opt$margin)
	keep <- which(ct$Contrast > opt$minContrast)
	if (length(keep) == 0) {
		msg("   %-28s 0/%d instances passed the contrast test, skipped", label, length(hits))
		return(list(name = label, status = "no_footprint", n = 0))
	}
	hits <- hits[keep]
	ct   <- ct[keep, , drop = FALSE]

	if (opt$maxInstances > 0 && length(hits) > opt$maxInstances) {
		top  <- order(ct$Contrast, decreasing = TRUE)[seq_len(opt$maxInstances)]
		hits <- hits[top]; ct <- ct[top, , drop = FALSE]
	}

	mcols(hits)$name <- sprintf("%s:%d", label, seq_along(hits))
	selBed <- data.frame(chr = as.character(seqnames(hits)),
						 start = start(hits) - 1, end = end(hits),
						 name = mcols(hits)$name, score = round(ct$Contrast, 4),
						 strand = as.character(strand(hits)),
						 Anchor = round(ct$Anchor, 5), Flank = round(ct$Flank, 5))
	## Note: this table is written once, as CnR.sorted.bed inside the motif's own
	## folder (same columns, ordered as plotted). The unfiltered scan for every
	## motif is in <outPrefix>.1.motifScan.bed.

	## Sequence and signal windows
	seqGr <- hits
	ranges(seqGr) <- IRanges(start = start(hits) - opt$seqMargin, end = end(hits) + opt$seqMargin)
	seqWin <- extractSeq(seqGr, faLocal)

	sigGr <- hits
	ranges(sigGr) <- IRanges(start = start(hits) - opt$exoMargin, end = end(hits) + opt$exoMargin)
	prof <- extractStranded1bp(sigGr, src.bwPlus, src.bwMinus, flipByStrand = TRUE)

	## Order instances by total signal, strongest first
	ord <- order(rowSums(cbind(prof$plus, -prof$minus)), decreasing = TRUE)
	seqWin <- seqWin[ord]
	fwd <- prof$plus[ord, , drop = FALSE]
	rev_ <- prof$minus[ord, , drop = FALSE]

	## Rank-prefixed so the most significant Homer motif sorts first. The rank is
	## zero-padded, otherwise motif10 would sort ahead of motif2.
	desDir <- sprintf("%s.3.motif%02d.%s", outPrefix, m$rank, label)
	dir.create(desDir, recursive = TRUE, showWarnings = FALSE)
	vizPrefix <- file.path(desDir, "CnR")

	writeXStringSet(seqWin, sprintf("%s.sorted.fa", vizPrefix))
	writeBed(selBed[ord, ], sprintf("%s.sorted.bed", vizPrefix))

	## Titles carry Homer's prediction for this motif
	homerTag <- sprintf("%s  [Homer motif%d: %s", m$name, m$rank, m$consensus)
	if (!is.na(m$pvalue)) homerTag <- sprintf("%s, p=%.1e", homerTag, m$pvalue)
	if (!is.na(m$target)) homerTag <- sprintf("%s, T=%.1f%%", homerTag, m$target)
	homerTag <- paste0(homerTag, "]")
	mainTitle <- sprintf("%s\n%s", dataName, homerTag)

	## Average profile alone
	renderFigure(sprintf("%s.viz.avg.png", vizPrefix), width = 5, height = 3.8,
				 density = opt$density, keepPdf = isTRUE(opt$pdf),
				 expr = {
					 par(las = 1, cex = 0.72)
					 visualizeAvgProfile(fwd, rev_, mainTitle, motifLen, yMaxAvg = yMaxAvg)
				 })

	## Sequence logos from the retained instances: motif locus only, and the
	## wider window used for the nucleotide heatmap. Both kept as PNG and PDF.
	motifSeq <- subseq(seqWin, opt$seqMargin + 1, opt$seqMargin + motifLen)
	renderFigure(sprintf("%s.logo.short.png", vizPrefix),
				 width = max(3.2, motifLen * 0.34), height = 2.4,
				 density = opt$density, keepPdf = TRUE,
				 expr = drawSeqLogo(motifSeq,
						sprintf("%s  motif only (%d bp, N=%d)", label, motifLen, length(hits))))
	renderFigure(sprintf("%s.logo.long.png", vizPrefix),
				 width = max(4.5, width(seqWin)[1] * 0.16), height = 2.4,
				 density = opt$density, keepPdf = TRUE,
				 expr = drawSeqLogo(seqWin,
						sprintf("%s  motif +/- %d bp (%d bp, N=%d)", label, opt$seqMargin,
								width(seqWin)[1], length(hits))))

	## Compact, title-less copy of the short logo for the combined figure's
	## upper-left panel (the standalone files above keep their titles).
	## Rendered at roughly the panel's own aspect so no vertical space is wasted,
	## then image_trim()ed to strip seqLogo's white border -- between them these
	## let the logo fill the panel instead of floating in the middle of it.
	logoRaster <- NULL
	embedPng <- tempfile(fileext = ".png")
	tryCatch({
		embedW <- max(3.2, motifLen * 0.34)
		renderFigure(embedPng, width = embedW, height = embedW * logoPanelAspect,
					 density = opt$density, keepPdf = FALSE,
					 expr = drawSeqLogo(motifSeq, "", fontsize = 8))
		logoRaster <- as.raster(magick::image_trim(magick::image_read(embedPng)))
	}, error = function(e) msg("   %-28s logo panel skipped (%s)", label, conditionMessage(e)))
	unlink(embedPng)

	## Combined three-panel figure
	fullTitle <- sprintf("%s\n(N=%d)", mainTitle, length(hits))
	renderFigure(sprintf("%s.viz.png", vizPrefix), width = 6, height = 8,
				 density = opt$density, keepPdf = isTRUE(opt$pdf),
				 expr = visualizeInstances(seqWin, fwd, rev_, title = fullTitle,
										   motifLen = motifLen, yMaxAvg = yMaxAvg,
										   yMaxHeatmap = yMaxHeatmap,
										   cex = 1, mode = opt$combineMode,
										   logoRaster = logoRaster))

	msg("   %-28s %d/%d instances kept -> %s", label, length(hits), length(ok), desDir)
	list(name = label, status = "ok", n = length(hits))
}

safeProcess <- function(m) {
	tryCatch(processMotif(m), error = function(e)
		list(name = m$name, status = paste("error:", conditionMessage(e)), n = 0))
}

nCore <- max(1L, as.integer(opt$parallel))
results <- NULL
if (nCore > 1 && .Platform$OS.type == "unix") {
	nCore <- min(nCore, detectCores())
	vmsg("Processing %d motif(s) on %d core(s)", length(selected), nCore)
	results <- suppressWarnings(mclapply(selected, safeProcess, mc.cores = nCore))
	## A forked worker that dies returns a try-error or NULL rather than our list;
	## fall back to serial so a fork problem never silently loses a motif.
	bad <- !vapply(results, function(r) is.list(r) && !is.null(r$status), logical(1))
	if (any(bad)) {
		msg("Warning: %d parallel worker(s) failed to return; re-running serially", sum(bad))
		results <- NULL
	}
}
if (is.null(results)) results <- lapply(selected, safeProcess)

## ---------------------------------------------------------------- summary ---
msg("")
summary.tab <- do.call(rbind, lapply(results, function(r)
	data.frame(Motif = r$name, Status = r$status, Instances = r$n, stringsAsFactors = FALSE)))
write.table(summary.tab, sprintf("%s.5.summary.txt", outPrefix),
			sep = "\t", quote = FALSE, row.names = FALSE)

failed <- grep("^error:", summary.tab$Status)
drawn  <- sum(summary.tab$Status == "ok")

msg("Summary:")
for (i in seq_len(nrow(summary.tab)))
	msg("  %-28s %-14s N=%d", summary.tab$Motif[i], summary.tab$Status[i], summary.tab$Instances[i])
msg("")

if (length(failed) > 0) {
	msg("Error: %d motif(s) failed", length(failed))
	quit(save = "no", status = 1)
}
if (drawn == 0) finishEmpty("No motif produced a footprint figure")

writeLines("", sprintf("%s.4.complete", outPrefix))
msg("Footprint analysis complete: %d motif(s) visualized", drawn)
