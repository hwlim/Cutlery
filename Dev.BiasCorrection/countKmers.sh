#!/usr/bin/env bash

trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT
source $COMMON_LIB_BASE/commonBash.sh

##############################################################
# Extract k-mer sequences around 5'-end of aligned reads
# and calculate frequency of unique k-mers
# If mask bed file is given, reads within the mask regions are considred only.


function printUsage {
	echo -e "Usage: `basename $0` (options) <input file>
Input:
	- BAM file or
	- BED file (gzipped) or
	- \"stdin\" to use /dev/stdin
Options:
	-l <k_up>: Upstream k-mer length, default=0
	-r <k_down>: Downstream k-mer lengt includg 5'-end, default=0
	-c <chrRegex>: chromosomes to consider in regular expression format, default=^chr[0-9XY]+$
	-g <genome>: genome (for Homer) or genome fasta file, required
	-s <chromSize>: chromosome size file, required
	-v : Verbose mode
Output:
	Two column output (k-mer and count) in STDOUT
	- Sorted by k-mers
	- Possibly some k-mer missing, i.e. meaning no zero count" >&2
}


###################################
## option and input file handling
k_up=0
k_down=0
chrRegex="^chr[0-9XY]+$"
genome=NULL
chromSize=NULL
verbose=FALSE
while getopts ":l:r:c:g:s:v" opt; do
	case $opt in
		l)
			k_up=$OPTARG
			;;
		r)
			k_down=$OPTARG
			;;
		c)
			chrRegex=$OPTARG
			;;
		g)
			genome=$OPTARG
			;;
		s)
			chromSize=$OPTARG
			;;
		v)
			verbose=TRUE
			;;
		\?)
			echo "Invalid options: -$OPTARG" >&2
			printUsage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			printUsage
			exit 1
			;;
	esac
done


if [ "$genome" = "NULL" ];then
	echo "Error: genome must be specified (-g)" >&2
	exit 1
fi

if [ "$chromSize" = "NULL" ];then
	echo "Error: chromosome size file must gbe specified (-s)" >&2
	exit 1
else
	assertFileExist $chromSize
fi

if [ $k_up -eq 0 ] && [ $k_down -eq 0 ];then
	echo -e "Error: Total k-mer length must be > 0" >&2
	exit 1
fi
kmerLen=`echo -e "${kmer_up}\t${kmer_dn}" | gawk '{ print $1 + $2 }'`



shift $((OPTIND-1))
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


src=$1
if [ "${src}" = "stdin" ];then
	src=/dev/stdin
else
	assertFileExist ${src}
fi



if [ "${verbose}" = "TRUE" ];then
	echo -e "==============================================">&2
	echo -e "Count K-mers at 5'-end" >&2
	echo -e "  Input = $src" >&2
	echo -e "  k-mer = $kmer_up / $kmer_dn" >&2
	echo -e "  genome = ${genome}" >&2
	echo -e "  chrRegex = ${chrRegex}" >&2
fi


## Input type classification
if [ "$src" = "/dev/stdin" ];then
	cmdCat=cat
else
	srcFile=`basename $src`
	ext=${srcFile##*.}
	if [ "${ext}" == "bam" ];then
		cmdCat="bamToBed -i"
	elif [ "${ext}" == "gz" ];then
		cmdCat=zcat
	else
		echo -e "Error: Invalid input type with an extension \"$ext\"" >&2
		exit 1
	fi
fi

$cmdCat $src \
	| gawk '$1 ~ /'$chrRegex'/ { if($6=="+"){ c=$2 }else{ c=$3 } if(c!=0) printf "%s\t%d\t%d\t%s\t0\t%s\n", $1,c,c,$4,$6 }' \
	| slopBed -i stdin -g $chromSize -s -l ${kmer_up} -r ${kmer_dn} \
	| gawk '$3-$2=='$kmer'' \
	| bedtools getfasta -fi $genomeFa -bed stdin -fo /dev/stdout -s -tab -name \
	| cut -f 2 \
	| tr '[a-z]' '[A-Z]' \
	| gawk '{
			if( $1 in cntDic ){
				cntDic[$1] = cntDic[$i] + 1
			}else{
				cntDic[$1] = 1
			}
		}
		END{
			for( seq in cntDic ) printf "%s\t%d\n", seq, cntDic[seq]
		}' \
	| sort -k1,1 
