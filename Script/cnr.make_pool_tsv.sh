#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh

tsv=$1

assertFileExist $tsv

grep -v ^# $tsv | head -n 1 | cut -f 1-7
grep -v ^# $tsv | grep -v ^$ | tail -n +2 | cut -f 3,7 | sort -k1,1 | uniq |  gawk '{ printf "%s\t%s\t%s\tNULL\tNULL\tNULL\t%s\n", $1, $1, $1, $2 }' 

