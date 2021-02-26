# Cutlery: automated pipeline and utilities for CUT&RUN data analysis

Near-future plan
- Remove "Pool" rules
Once bam files are pooled, all the downstream analysis is the same
Therefore, it is more logical and simple to have one script for bam file pooling utilizig sample.tsv


## 0. Prerequisite

### 0.1. LimLabBase

This pipeline heavily relies on LimLabBase. Therefore, setup github.com/LimLabBase first.

```bash
git clone https://github.com/LimLabBase
```
And add the subfolders of the LimLabBase in the PATH of .bash_profile or .bashrc

### 0.2. Other software packages
- Homer
- cutadapt
- STAR (> v2.7.4)
- bedtools
- samtools
- bedGraphToBigWig (UCSC utilities)
- bwtool
- gawk
- sed
- R
- 
All available in the HPC/BMI


## 1. Setting Cutlery for CUT&RUN data analysis

### 1.1. Git clone:
```bash
git clone https://github.com/hwlim/Cutlery
```
Currently, **LimLabBase** also contains a folder, Cutlery. But **Cutlery** can be clonned to any location as far as you properly declare the environment variable below.

### 1.2. Declare environment for Cutlery in .bash_profile (or .bashrc)

```bash
# For Cutlery location
# This environment variable is used in many Cutlery scripts to locate basic resource
# Replace <path_to_Cutlery> with as needed
export CUTLERY=<path_to_Cutlery>

# PATH for Cutlery scripts and executables
export PATH=${PATH}:${CUTLERY}
```

## 2. Analysis of CUT&RUN data

### 2.1. Initialize analysis folder

Create a folder for analysis workspace

```bash
mkdir CnR_Analysis
cd CnR_Analysis
```

Initialize to create template files

```bash
cnr.init.sh
```
which will create three files:
- sample.tsv: sample information meta sheet
- Snakefile: analysis parameters and environments
- 0.submit.snakemake.sh: script to run the analysis


### 2.2. sample.tsv

Tab-separated sample information file with following 7 columns:
- **Id**: Unique sample ID. This is used for the output files of adapter trimming. 
- **Name**: Sample name. This becomes the output folder name for each sample
- **Group**: Sample group. No specific use in the process of each sample but will be used for replicate-pooling and group-wise analysis later
- **Fq1 & 2**: Name of fastq files, Read1 & 2 (Paired-end sequencing is assumed always)
- **Ctrl**: Name of control sample from the "Name" column for peak calling. Use "NULL" for peak calling without a control or ATAC-seq data analysis.
- **PeakMode**: Peak mode among {factor,histone,NULL}. If "NULL", no peak calling is performed, e.g. for IgG samples.

Note:
- No white space is allowed in the file except for the column separator, tab.
- **Id** and **Name** columns must contain unique values not overlapping with other samples

Example:

|Id|Name|Group|Fq1|Fq2|Ctrl|PeakMode|
|--|----|-----|---|---|----|--------|
|sample1|hPSC_Foxa1_rep1|hPSC_Foxa1|sample1_R1.fq.gz|sample1_R2.fq.gz|hPSC_IgG|factor|
|sample2|hPSC_Foxa1_rep2|hPSC_Foxa1|sample2_R1.fq.gz|sample2_R2.fq.gz|hPSC_IgG|factor|
|sample3|hPSC_IgG|hPSC_IgG|sample3_R1.fq.gz|sample3_R2.fq.gz|NULL|NULL|
|sample4|hPSC_H3K27me3_rep1|hPSC_H3K27me3|sample4_R1.fq.gz|sample4_R2.fq.gz|hPSC_IgG|histone|
|sample5|hPSC_H3K27me3_rep2|hPSC_H3K27me3|sample5_R1.fq.gz|sample5_R2.fq.gz|hPSC_IgG|histone|


### 2.3. Snakefile
- **Snakefile** declares various parameters, environments, files, and outputs.
- Check the comments within the **Snakefile** for further information and revise as needed.
- Basically, it is based on python. Mind the grammar.

## 3. Run
Dry run first to see if everything is correctly defined.
Snakemake is available in python3/3.6.3 in CCHMC/HPC.
```bash
module load python3/3.6.3
snakemake -np
```
If there's an error, find the source of error and correct them such as incorrect parameters, folders, or missing inputs.
If no error, then submit a Snakemake job for actual analysis.
```bash
./0.submet.snakemake.sh
```

After submitting, check if jobs are automatically created and submitted by Snakemake
```bash
bjobs
```
During the first few seconds, there will be only one job, a master job that create/submit each analysis task.
But eventually a few seconds later, children bjos will appear and run as their dependency is satisfied.
There is no automatic notification. Check back with **bjobs** command if they are still running.

## 4. Check results
Once the job is finished or encounter an error, two log files will be created.
(These file names are specified in 0.submit.snakemake.sh)
- bsub.err
- bsub.out

Try
```bash
tail bsub.err
```
If you see something like **"(100%) done"**, everything is OK and complete.
If you see something less than 100%, it means there was an error in one of the job.
Check the log of individual jobs under the folder **logs** to see the error message.
```baseh
grep -i error logs/*err
```
