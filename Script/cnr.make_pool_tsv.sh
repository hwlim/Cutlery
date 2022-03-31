#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

tsv=$1

assertFileExist $tsv

head -n 1 $tsv
tail -n +2 $tsv | cut -f 3,7 | sort -k1,1 | uniq |  gawk '{ printf "%s\t%s\t%s\tNULL\tNULL\tNULL\t%s\n", $1, $1, $1, $2 }' 

