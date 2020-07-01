#/usr/bin/env bash


########################################################################################
# Calculate average base frequency in background reads for both of ChIP-exo and seq

trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPIR}/__temp__.$$.*; fi' EXIT
source commonBash.sh


desDirSuffix=BaseFrequency

checkBaseFreq(){
	local src=$1
	local prefix=$2

	local desFreq=${desDir}/${prefix}.freq.txt
	local desPlotPrefix=${desDir}/${prefix}.freq.line
	local desPlot=${desPlotPrefix}.png

	local name=`dirname $workDir`
	name=`basename $name`

	echo -e "Drawing base frequency" >&2
	echo -e "  src = $src ">&2
	echo -e "  genome = $genome" >&2
	echo -e "  chrRx = ${chrRx}" >&2
	echo -e "  desFreq = $desFreq" >&2
	echo -e "  desPlot = ${desPlot}" >&2

	if [ -f $desFreq ];then
		echo -e "  - $desFreq already exists, pass" >&2
	else
		echo -e "  - Checking nucleotide base frequency: $src" >&2
		zcat $src \
			| gawk '/^'$chrRx'/' \
			| checkBaseFreq.sh -g $genome -l 20,20 -o ${desFreq} stdin
	fi

	if [ ! -f $desPlot ] || [ $desPlot -ot $desFreq ];then
		~/Research/ExoTools/7.ExoKmer.autosome.TSVrd/drawFreqLinePlot.r -x Offset -y Frequency -t "${name}.${prefix}" -l "0.0,0.5" -s 8,4 \
			-c firebrick1,royalblue,gold,limegreen -o $desPlotPrefix $desFreq
	else
		echo -e "  - $desPlot already exits, pass" >&2
	fi
	echo -e "" >&2
}

startDir=`pwd`
for workDir in `ls -d */Background/`
do
	echo -e "================================================" >&2
	echo -e "Entering $workDir" >&2
	cd $workDir
	if [ -f ./config.sh ];then
		echo -e "  - Loading config.sh" >&2
		source ./config.sh
	else
		continue
	fi

	desDir=${desDirSuffix}
	mkdir -p $desDir

	if [ -f ${exo} ];then
		checkBaseFreq ${exo} exo
	else
		echo -e "Warning: $exo does not exists, pass" >&2
	fi

	if [ -f ${chip} ];then
		checkBaseFreq ${chip} seq
	else
		echo -e "Warning: $chip does not exists, pass" >&2
	fi

	cd $startDir
	echo "" >&2
done

