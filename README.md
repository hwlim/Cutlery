# Cutlery: automated pipeline and utilities for CUT&RUN data analysis


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

## 2. Analysis of CUT&RUN data

### 2.1 sample.tsv file: sample data sheet

|Id|Name|Group|Fq1|Fq2|Ctrl|PeakMode|
|--|----|-----|---|---|----|--------|
|sample1|hPSC_Foxa1_rep1|hPSC_Foxa1|sample1_R1.fq.gz|sample1_R2.fq.gz|hPSC_IgG|factor|
|sample2|hPSC_Foxa1_rep2|hPSC_Foxa1|sample2_R1.fq.gz|sample2_R2.fq.gz|hPSC_IgG|factor|
|sample3|hPSC_IgG|hPSC_IgG|sample3_R1.fq.gz|sample3_R2.fq.gz|hPSC_IgG|NULL|
|sample4|hPSC_H3K27me3_rep1|hPSC_H3K27me3|sample4_R1.fq.gz|sample4_R2.fq.gz|hPSC_IgG|histone|
|sample5|hPSC_H3K27me3_rep2|hPSC_H3K27me3|sample5_R1.fq.gz|sample5_R2.fq.gz|hPSC_IgG|histone|
