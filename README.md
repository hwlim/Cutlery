# Cutlery: automated pipeline and utilities for CUT&RUN data analysis

## 0. Prerequisite

1. LimLabBase

This pipeline heavily relies on LimLabBase. Therefore, setup github.com/LimLabBase first.

```bash
git clone https://github.com/LimLabBase
```

And add the subfolders of the LimLabBase in the PATH of .bash_profile or .bashrc


## 1. Setting Cutlery for CUT&RUN data analysis

1. Git clone:
```bash
git clone https://github.com/hwlim/Cutlery
```

2. Declare environment for Cutlery in .bash_profile (or .bashrc)

```bash
# For Cutlery location
# This environment variable is used in many Cutlery scripts to locate basic resource
export CUTLERY=${HOME}/bin/Cutlery

# PATH for Cutlery scripts and executables
export PATH=${PATH}:${CUTLERY}
```

## 2. Analysis of CUT&RUN data: Each sample

### 2.0 Initialize analysis folder

Create a folder for analysis workspace

```bash
mkdir CnR.20201122
cd CnR.20201122
```

Initialize to create template files

```bash
cnr.init.sh
```
which will create three files:
- sample.tsv
- Snakefile
- 0.submit.snakemake.sh


### 2.1 sample.tsv file: sample data sheet with 6 columns:

- **Id**: Unique sample ID. This is used for the output files of adapter trimming. 
- **Name**: Sample name. This becomes the output folder name for each sample
- **Group**: Sample group. No specific use in the process of each sample but will be used for replicate-pooling and group-wise analysis later
- **Fq1 & 2**: Name of fastq files, Read1 & 2 (Paired-end sequencing is assumed always)
- **Ctrl**: Name of control sample from the "Name" column for peak calling. Use "NULL" for peak calling without a control.
- **PeakMode**: Peak mode among {factor,histone,NULL}. If "NULL", no peak calling is performed, e.g. for IgG samples.

Note: no white space is allowed in the file except for the column separator, tab.

Example:

|Id|Name|Group|Fq1|Fq2|Ctrl|PeakMode|
|--|----|-----|---|---|----|--------|
|sample1|hPSC_Foxa1_rep1|hPSC_Foxa1|sample1_R1.fq.gz|sample1_R2.fq.gz|hPSC_IgG|factor|
|sample2|hPSC_Foxa1_rep2|hPSC_Foxa1|sample2_R1.fq.gz|sample2_R2.fq.gz|hPSC_IgG|factor|
|sample3|hPSC_IgG|hPSC_IgG|sample3_R1.fq.gz|sample3_R2.fq.gz|hPSC_IgG|NULL|
|sample4|hPSC_H3K27me3_rep1|hPSC_H3K27me3|sample4_R1.fq.gz|sample4_R2.fq.gz|hPSC_IgG|histone|
|sample5|hPSC_H3K27me3_rep2|hPSC_H3K27me3|sample5_R1.fq.gz|sample5_R2.fq.gz|hPSC_IgG|histone|


