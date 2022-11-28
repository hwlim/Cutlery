#!/usr/bin/env Rscript

## NOTE: K-mer correction parameters are not currently validated

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('yaml', quiet=TRUE))


source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))

# command line option handling
option_list <- list(
	#make_option(c("-o","--outPrefix"), default=NULL, help="Prefix of output files, default=< anchor file name except for path and extension")
)
parser <- OptionParser(usage = "%prog [config_cnr.yml]", option_list=option_list,
	description="Validate given CUT&RUN config file (YAML)
Input:
	- CUT&RUN config file in YAML format, default=./config_cnr.yml" )

#################################
# Option handling
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) < 1){
	src="./config_cnr.yml"
}else{
	src=arguments$args[1]
}

opt=arguments$opt
#outPrefix=opt$outPrefix
#margin=opt$margin
#pseudo=opt$pseudo
#offset=opt$offset

assertFileExist(src)
config = yaml.load_file(src)


write(sprintf("Validating Cutlery config file: %s", src), stderr())
write(sprintf("Checking essential files"), stderr())
essentialFiles=c(
	"src_sampleInfo",
	"genomeFa",
	"chrom_size"
)
for( name in essentialFiles ){
	if(name %in% names(config)){
		assertFileExist(config[[name]])
	}else{
		stop(sprintf("<%s> must be  specified in %s", name, src))
	}
}

write(sprintf("Checking optional files"), stderr())
optionalFiles=c( "peak_mask", "cluster_yml" )
for( name in optionalFiles ){
	if(name %in% names(config[[name]])){
		assertFileExist(name)
	}
}

write(sprintf("Checking essential folders"), stderr())
essentialDirs=c( "star_index" )
for( name in essentialDirs ){
	if(name %in% names(config)){
		assertDirExist(config[[name]])
	}else{
		stop(sprintf("<%s> must be  specified in %s", name, src))
	}
}

write(sprintf("Checking essential string parameters"), stderr())
essentialChar=c(
	"genome",
	"adapter",
	"chrRegexAll",
	"chrRegexTarget",
	"star_option"
)
for( name in essentialChar ){
	if(name %in% names(config)){
		if(!is.character(config[[name]])) stop(sprintf("<%s> must be string in %s", name, src))
	}else{
		stop(sprintf("<%s> must be specified in %s", name, src))
	}
}

write(sprintf("Checking essential numeric parameters"), stderr())
essentialNumeric=c(
	"trim_minLen",
	"trim_minQual"
)
for( name in essentialNumeric ){
	if(name %in% names(config)){
		if(!is.numeric(config[[name]])) stop(sprintf("<%s> must be numeric in %s", name, src))
	}else{
		stop(sprintf("<%s> must be specified in %s", name, src))
	}
}

src_sampleInfo  = config[["src_sampleInfo"]]
write(sprintf("Checking sample information: %s", src_sampleInfo), stderr())
samples = read.delim(src_sampleInfo, header=TRUE, stringsAsFactors=FALSE, comment.char="#")
hdr = c("Id","Name","Group","Fq1","Fq2","Ctrl","PeakMode")
if( any(colnames(samples) != hdr) ) stop(sprintf("Invalid sample table format in %s. Must include columns: %s", src_sampleInfo, paste(hdr, collapse=" / ")))
if( nrow(samples) != length(unique(samples$Id)) ) stop( "Id column must contains unique values" )
if( nrow(samples) != length(unique(samples$Name)) ) stop( "Name column must contains unique values" )
tmp = setdiff( samples$Ctrl, c( samples$Name, "NULL" ) )
if( length(tmp) > 0 ) stop( sprintf("Unknown name in Ctrl column (must match with Name column): %s", paste(tmp, collapes=" / ")) )
tmp = setdiff( samples$PeakMode, c( "factor", "histone", "NULL" ) )
if( length(tmp) > 0 ) stop( sprintf("Invalid values in PeakMode column (must be either factor / histone / NULL): %s", paste(tmp, collapes=" / ")) )

write("Done!", stderr())

#kmer_genome     = config["kmer_genome"]
#kmer_pseudo     = config["kmer_pseudo"]

