#!/usr/bin/env python3

'''
    Written by:
        Christopher Ahn
        Hee Woong Lim

    Script to convert coordinate-sorted BAM file (with index) into a bed file of fragment or read.
    - fragment: interconnecting read1 & read2 from paired-ends
    - read: treat read1 & read2 (if any) separately
'''

#import required modules and check if user has pysam installed
try:
    import pysam
except ImportError:
    print("Error: This script requires the module \"pysam\" to run.")
    exit()
import argparse
import os
import re

#parse command line arguments
options = argparse.ArgumentParser(description="Converts a coordinate-sorted bam file to a fragment bed format. 5th column of the output represents the map quality.", usage="python sortedBamToFrag.py (options) [bam]")
options.add_argument('-f', '--flags_include', default='0x2',
                        help='SAM flag to include; input can be either in decimal or hexadecimal format. Default = 0x2. Use \'-f NULL\' to include all SAM flags (this means no filtering; -F should be set to NULL as well in this case). User can enter multiple flags by entering their sum; ex. if the user wants to include flags 2 and 64, type \"-f 66\" or \"-f 0x42\" without the quotation marks.')
options.add_argument('-F', '--flags_exclude', default='0x400',
                        help='SAM flag to exclude; input can be either in decimal or hexadecimal format. Default = 0x400. Use \'-F NULL\' along with \'-f NULL\' to include all SAM flags (no exclusion). User can enter multiple flags by entering their sum; ex. if the user wants to exclude flags 512 and 1024, type \"-F 1536\" or \"-f 0x600\" without the quotation marks.')
options.add_argument('-c', '--chr_include', default='.',
                        help='Regular expression of chromosomes to select. Default = . (all). e.g.) ^chr[0-9XY]+$|^chrM$ : regular/sex/chrM, ^chr[0-9XY]+$ : autosomal and sex chromosomes only.')
group = options.add_mutually_exclusive_group()
group.add_argument('-p', action='store_true',
                        help='Print a pair of reads in one line; same functionality as bamToBed -bedpe. Add the -p flag to use. Cannot be used with the -r flag. Using this option will not create fragments.')
group.add_argument('-r', action='store_true',
                        help='Convert bam file to bed format; same functionality as bedtools bamToBed with default settings, i.e. assuming single-end. Add the -r flag to use. Cannot be used with the -p flag. Using this option will not create fragments.')
options.add_argument('bam_file',
                        help='Required; Coordinate sorted bam file. Index file needs to be in the same location as the bam file.')
args = options.parse_args()

#function to extract read length from CIGAR string (column 6 in bam file)
def getReadLen(cig):
    #define read length variable
    readLen = 0
    
    #Include both deletions and matches in the total read length
    #need to loop through entire CIGAR string as there can be indels between matches.
    for tup in cig:
        #matches are labeled as "0", and deletions are labeled as "2" by pysam
        if (tup[0] == 0) or (tup[0] == 2):
            readLen += tup[1]
    
    return readLen

#function that performs bamToBed functionality, i.e. single-end mode
def bamToBed(bamFile, chrPattern):

    if args.flags_include != "NULL":
        #parse command line flags and convert hexadecimal inputs accordingly
        try:
            int(args.flags_include)
            inFlags = int(args.flags_include)
        except ValueError:
            inFlags = int(args.flags_include, 16)
    
    else:
        inFlags = 0

    if args.flags_exclude != "NULL":
        try:
            int(args.flags_exclude)
            exFlags = int(args.flags_exclude)
        except ValueError:
            exFlags = int(args.flags_exclude, 16)
    else:
        exFlags = 4096

    #Read the bam file line-by-line
    for read in bamFile.fetch():

        #Check flags and chromosome and skip reads that don't fit user specified criteria
        if (int(read.flag) & inFlags == inFlags) and (int(read.flag) & exFlags != exFlags) and (chrPattern.match(bamFile.get_reference_name(read.reference_id))) :

            #get chromosome name, as chromosomes are encoded as integers in sam/bam file
            chromName = bamFile.get_reference_name(read.reference_id)
            
            #get read length, start position, end position from current line
            readLen = getReadLen(read.cigartuples)
            start = int(read.pos)
            end = start + readLen

            #Get strand direction from column 9 in bam file
            if int(read.tlen) > 0:
                strand = "+"
            else:
                strand = "-"

            #check if the current line is the first segment or last segment of the read
            if int(read.flag) & 64 == 64 :
                tail = "/1"
            elif int(read.flag) & 128 == 128 :
                tail = "/2"

            #get read name and add suffix according to order of segment in template
            readName = read.qname + tail

            #print read
            printLine = [chromName, start, end, readName, read.mapq, strand]
            print(*printLine, sep="\t")
    return


#function that converts coordinate-sorted bam file to a fragment.bed file, i.e. paired-end mode
def bamToFrag(bamFile, chrPattern):
    #create empty dictionary
    d = {}

    if args.flags_include != "NULL":
        #parse command line flags and convert hexadecimal inputs accordingly
        try:
            int(args.flags_include)
            inFlags = int(args.flags_include)
        except ValueError:
            inFlags = int(args.flags_include, 16)
    
    else:
        inFlags = 0

    if args.flags_exclude != "NULL":
        try:
            int(args.flags_exclude)
            exFlags = int(args.flags_exclude)
        except ValueError:
            exFlags = int(args.flags_exclude, 16)
    else:
        exFlags = 4096

    #iterate through each line in bam file
    for read in bamFile.fetch():

        #Check flags and chromosome and skip reads that don't fit user specified criteria
        if (int(read.flag) & inFlags == inFlags) and (int(read.flag) & exFlags != exFlags) and (chrPattern.match(bamFile.get_reference_name(read.reference_id))) :

            #get read name from current line
            readName = read.qname

            #if read name exists in dictionary, pop key and get value
            if readName in d:
                firstEncounteredRead = d.pop(readName, None)

                #assign relevant variables to the values after popping the key
                firstPOS = firstEncounteredRead[0]
                firstSTRAND = firstEncounteredRead[1]
                
                #obtain chromosome ID from current line
                chromName = bamFile.get_reference_name(read.reference_id)

                #define start and end according to strand direction
                #if the first read is a "+" strand, the end position would be read_length + start_position of current line in the bam file
                if firstSTRAND == "+":
                    readLen = getReadLen(read.cigartuples)
                    end = int(read.pos) + readLen
                    start = firstPOS

                #build fragment as if it were a "+" strand if the first read seen was a "-"
                else:
                    end = firstPOS
                    start = int(read.pos)
                    firstSTRAND = "+"

                #print fragment info to STDOUT
                readInfo = [chromName, start, end, readName, read.mapq, firstSTRAND]
                print(*readInfo, sep="\t")
            

            #if read name doesn't exist in dictionary, add relevant info to the dictionary
            else:
                #get read length of current line
                readLen = getReadLen(read.cigartuples)

                #get strand direction using TLEN column in bam file
                if int(read.tlen) > 0:
                    strand = "+"
                    POS = int(read.pos)
                
                #if the strand is "-", the end position would be the starting position + read length of the current line in the bam file
                else:
                    strand = "-"
                    POS = int(read.pos) + readLen

                #create key and value for new unseen read and add to dictionary
                d[readName] = [POS, strand]
    return

def bedpe(bamFile, chrPattern):
    #iterate through each line in bam file
    d = {}
    
    if args.flags_include != "NULL":
        #parse command line flags and convert hexadecimal inputs accordingly
        try:
            int(args.flags_include)
            inFlags = int(args.flags_include)
        except ValueError:
            inFlags = int(args.flags_include, 16)
    
    else:
        inFlags = 0

    if args.flags_exclude != "NULL":
        try:
            int(args.flags_exclude)
            exFlags = int(args.flags_exclude)
        except ValueError:
            exFlags = int(args.flags_exclude, 16)
    else:
        exFlags = 4096
        
    for read in bamFile.fetch():

        #Check flags and chromosome and skip reads that don't fit user specified criteria
        if (int(read.flag) & inFlags == inFlags) and (int(read.flag) & exFlags != exFlags) and (chrPattern.match(bamFile.get_reference_name(read.reference_id))) :

            #get read name from current line
            readName = read.qname

            #if read name exists in dictionary, pop key and get value
            if readName in d:
                firstEncounteredRead = d.pop(readName, None)

                #assign relevant variables to the values after popping the key
                firstStart = firstEncounteredRead[0]
                firstEnd = firstEncounteredRead[1]
                firstSTRAND = firstEncounteredRead[2]
                firstChrom = firstEncounteredRead[3]
                
                #work on current read
                readLen = getReadLen(read.cigartuples)
                currentStart = int(read.pos)
                currentEnd = int(read.pos) + readLen
                currentChrom = bamFile.get_reference_name(read.reference_id)

                #define start and end according to strand direction
                #if the first read is a "+" strand, the end position would be read_length + start_position of current line in the bam file
                if firstSTRAND == "+":
                    currentSTRAND = "-"

                #build fragment as if it were a "+" strand if the first read seen was a "-"
                else:
                    currentSTRAND = "+"

                    
                #print fragment info to STDOUT
                readInfo = [firstChrom, firstStart, firstEnd, currentChrom, currentStart, currentEnd, readName, read.mapq, firstSTRAND, currentSTRAND]
                print(*readInfo, sep="\t")


            #if read name doesn't exist in dictionary, add relevant info to the dictionary
            else:
                #get read length of current line
                readLen = getReadLen(read.cigartuples)

                #get strand direction using TLEN column in bam file
                if int(read.tlen) > 0:
                    strand = "+"

                #if the strand is "-", the end position would be the starting position + read length of the current line in the bam file
                else:
                    strand = "-"
                
                start = int(read.pos)
                end = int(read.pos) + readLen
                chromName = bamFile.get_reference_name(read.reference_id)

                #create key and value for new unseen read and add to dictionary
                d[readName] = [start, end, strand, chromName]
    return

#main function
def main():
    #read input bam file
    bamFile = pysam.AlignmentFile(args.bam_file, "rb")
    chrPattern = re.compile(str(args.chr_include))

    if args.p:
        ##Print a pair of reads if user uses -p flag
        bedpe(bamFile, chrPattern)
    elif args.r:
        ##Run bamToBed if user uses -r flag
        bamToBed(bamFile, chrPattern)
    else:
        ##Run bamToFrag otherwise
        bamToFrag(bamFile, chrPattern)

if __name__ == "__main__":
    main()