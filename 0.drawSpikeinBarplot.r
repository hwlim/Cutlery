#!/usr/bin/env Rscript

library(ggplot2)
library(cowplot)
library(Seurat)
source(sprintf("%s/commonR.r", Sys.getenv("COMMON_LIB_BASE")))
source(sprintf("%s/basicR.r", Sys.getenv("COMMON_LIB_BASE")))

sampleInfo = read.delim("../../sample.tsv", header=TRUE, stringsAsFactors=FALSE)
sampleL = sampleInfo[,2]

spikeCntL = NULL

tmp.df = NULL
for( sample in sampleL ){
	# sample=sampleL[1]
	tmp = read.delim(sprintf("%s.spikeCnt.txt", sample), header=TRUE, row.names=1)

	if(is.null(tmp.df)){
		tmp.df = tmp[1:2, ,drop=F]
	}else{
		tmp.df = cbind(tmp.df, tmp[1:2,])
	}
}
colnames(tmp.df) = sampleL

#df = as.data.frame(t(tmp.df))
df = data.frame(Sample = sampleL, Target = as.numeric(tmp.df[1,]), Spike = as.numeric(tmp.df[2,]))
rownames(df) = NULL
df$Sample = factor(df$Sample, levels=sampleL)
df$TargetSpikeRatio = df$Target / df$Spike
df$SpikeFraction = df$Spike / (df$Target + df$Spike) * 100

if(FALSE){
df$Condition = sapply(sampleL, function(x) strsplit(x, "_")[[1]][3])
df$Condition = factor(df$Condition, levels=c("Control", "FOXA123KD"))
df$ChIP = sapply(sampleL, function(x) strsplit(x, "_")[[1]][2])
df$ChIP = factor(df$ChIP, levels=c("Input", "H3K9me3"))
df$Sample2 = sapply(sampleL, function(x) paste(strsplit(x, "_")[[1]][c(1,4)], collapse="_"))
df$Sample2 = factor(df$Sample2, levels=unique(df$Sample2))
}

g.Target = ggplot(data=df, aes(x=Sample, y=Target)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
	labs(title="Target", y="Read depth")


g.Spike = ggplot(data=df, aes(x=Sample, y=Spike)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
	labs(title="Spike", y="Read depth")

g.frac = ggplot(data=df, aes(x=Sample, y=SpikeFraction)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),) +
	ylab("Spikein Fraction (%)")

g.ratio = ggplot(data=df, aes(x=Sample, y=TargetSpikeRatio)) +
	geom_bar(stat="identity", fill="steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),) +
	ylab("Reads Ratio: Target / Spike")


#	g1 = plot_grid(g.Target, g.Spike, ncol=1)
g = plot_grid(g.Target, g.Spike, g.frac, g.ratio, ncol=2)

ggsave("TargetSpikeRatio.png", g, width=3+0.3*nrow(df), height=8, dpi=200)




